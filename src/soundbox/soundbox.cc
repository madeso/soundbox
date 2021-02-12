#include "soundbox/soundbox.h"

#include <cmath>
#include <cassert>
#include <unordered_map>
#include <functional>


namespace
{
    template<typename K, typename V, typename F>
    const auto& find_or_create(std::unordered_map<K, V>* map, const K& k, F&& f)
    {
        auto found = map->find(k);
        if (found != map->end())
        {
            return found->second;
        }
        else
        {
            auto where = map->emplace(std::make_pair(k, f()));
            return where.first->second;
        }
    }


    int cint(double d)
    {
        return static_cast<int>(d);
    }

    int cint(std::size_t i)
    {
        return static_cast<int>(i);
    }


    float d2f(double d)
    {
        return static_cast<float>(d);
    }


    auto generate_seed()
    {
        std::random_device device;
        return device();
    }


    float generate_random(Player::RandomEngine* random_engine)
    {
        auto distribution = std::uniform_real_distribution<float>{ 0.0f, 1.0f };
        return distribution(*random_engine);
    }

    i32 truncate(float f)
    {
        return static_cast<i32>(std::floor(f));
    }
}


namespace
{
    //--------------------------------------------------------------------------
    // Private methods
    //--------------------------------------------------------------------------


    float frac(float f)
    {
        return f - std::trunc(f);
    }


    // Oscillators
    float osc_sin(float value)
    {
        return d2f(std::sin(value * 6.283184));
    };


    float osc_saw(float value)
    {
        return 2 * frac(value) - 1;
    };


    float osc_square(float value)
    {
        return frac(value) < 0.5 ? 1.0f : -1.0f;
    };


    float osc_tri(float value)
    {
        auto v2 = frac(value) * 4;
        if (v2 < 2) return v2 - 1;
        return 3 - v2;
    };


    // Array of oscillator functions
    float mOscillators(int oscilator_type, float value)
    {
        switch (oscilator_type)
        {
        case 0:  return osc_sin(value);
        case 1:  return osc_square(value);
        case 2:  return osc_saw(value);
        case 3:  return osc_tri(value);
        default: assert(false); return 0.0f;
        }
    };


    float getnotefreq(float n)
    {
        // 174.61.. / 44100 = 0.003959503758 (F3)
        return 0.003959503758f * std::pow(2.0f, (n - 128) / 12.0f);
    };


    std::vector<i32> createNote(const Instrument& instr, int n, int rowLen, Player::RandomEngine* random_engine)
    {
        auto osc1 = [&](float f) { return mOscillators(instr.i[OSC1_WAVEFORM], f); };
        auto o1vol = instr.i[OSC1_VOL];
        auto o1xenv = instr.i[OSC1_XENV] / 32.0f;
        auto osc2 = [&](float f) { return mOscillators(instr.i[OSC2_WAVEFORM], f); };
        auto o2vol = instr.i[OSC2_VOL];
        auto o2xenv = instr.i[OSC2_XENV] / 32.0f;
        auto noiseVol = instr.i[NOISE_VOL];
        auto attack = instr.i[ENV_ATTACK] * instr.i[ENV_ATTACK] * 4;
        auto sustain = instr.i[ENV_SUSTAIN] * instr.i[ENV_SUSTAIN] * 4;
        auto release = instr.i[ENV_RELEASE] * instr.i[ENV_RELEASE] * 4;
        auto releaseInv = 1.0f / release;
        auto expDecay = -instr.i[ENV_EXP_DECAY] / 16.0f;
        auto arp = instr.i[ARP_CHORD];
        auto arpInterval = rowLen * cint(std::pow(2, (2 - instr.i[ARP_SPEED])));

        auto noteBuf = std::vector<i32>(attack + sustain + release);

        // Re-trig oscillators
        float c1 = 0, c2 = 0;

        // Local variables.
        int j, j2;
        float o1t=0.0f, o2t=0.0f;
        float e, rsample;

        // Generate one note (attack + sustain + release)
        for (j = 0, j2 = 0; j < attack + sustain + release; j++, j2++) {
            if (j2 >= 0) {
                // Switch arpeggio note.
                arp = (arp >> 8) | ((arp & 255) << 4);
                j2 -= arpInterval;

                // Calculate note frequencies for the oscillators
                o1t = getnotefreq(n + (arp & 15) + instr.i[OSC1_SEMI] - 1280.f);
                o2t = getnotefreq(n + (arp & 15) + instr.i[OSC2_SEMI] - 128.0f) * (1.0f + 0.0008f * instr.i[OSC2_DETUNE]);
            }

            // Envelope
            e = 1;
            if (j < attack) {
                e = static_cast<float>(j) / attack;
            }
            else if (j >= attack + sustain) {
                e = (j - attack - sustain) * releaseInv;
                e = (1 - e) * std::pow(3.0f, (expDecay * e));
            }

            // Oscillator 1
            c1 += o1t * std::pow(e, o1xenv);
            rsample = osc1(c1) * o1vol;

            // Oscillator 2
            c2 += o2t * std::pow(e, o2xenv);
            rsample += osc2(c2) * o2vol;

            // Noise oscillator
            if (noiseVol) {
                rsample += (2 * generate_random(random_engine) - 1) * noiseVol;
            }

            // Add to (mono) channel buffer
            noteBuf[j] = truncate(80 * rsample * e);
        }

        return noteBuf;
    }
}


//--------------------------------------------------------------------------
// Initialization
//--------------------------------------------------------------------------

Player::Player(const Song& song)
    : random_engine(generate_seed())
{
    // Define the song
    mSong = song;

    // Init iteration state variables
    mLastRow = song.endPattern;
    mCurrentCol = 0;

    // Prepare song info
    mNumWords =  song.rowLen * song.patternLen * (mLastRow + 1) * 2;

    // Create work buffer (initially cleared)
    mMixBuf.resize(mNumWords);
};


//--------------------------------------------------------------------------
// Public methods
//--------------------------------------------------------------------------

// Generate audio data for a single track
float Player::generate()
{
    // Local variables
    std::optional<int> n;
    int i, j, p, row, col,
        k, rowStartSample;

    float f, t;

    // Put performance critical items in local variables
    auto chnBuf = std::vector<i32>(mNumWords);
    auto& instr = mSong.songData[mCurrentCol];
    int rowLen = mSong.rowLen;
    int patternLen = mSong.patternLen;

    // Clear effect state
    float low = 0, band = 0, high;
    float lsample;
    float rsample;
    bool filterActive = false;

    // Clear note cache.
    auto noteCache = std::unordered_map<int, std::vector<i32>>();

    // Patterns
    for (p = 0; p <= mLastRow; ++p) {
        auto& cp = instr.p[p];

        // Pattern rows
        for (row = 0; row < patternLen; ++row) {
            // Execute effect command.
            auto cmdNo = cp ? instr.c[cp.value() - 1].f[row] : 0;
            if (cmdNo) {
                instr.i[cmdNo.value() - 1] = instr.c[cp.value() - 1].f[row + patternLen] || 0;

                // Clear the note cache since the instrument has changed.
                if (cmdNo < 17) {
                    noteCache.clear();
                }
            }

            // Put performance critical instrument properties in local variables
            auto oscLFO = [&](float f) { return mOscillators(instr.i[LFO_WAVEFORM], f); };
            auto lfoAmt = instr.i[LFO_AMT] / 512.0f;
            auto lfoFreq = static_cast<float>(std::pow(2.0f, (instr.i[LFO_FREQ] - 9))) / rowLen;
            auto fxLFO = instr.i[LFO_FX_FREQ];
            auto fxFilter = instr.i[FX_FILTER];
            auto fxFreq = instr.i[FX_FREQ] * 43.23529f * 3.141592f / 44100.0f;
            auto q = 1 - instr.i[FX_RESONANCE] / 255.0f;
            auto has_dist = instr.i[FX_DIST] != 0;
            auto dist = instr.i[FX_DIST] * 1e-5f;
            auto drive = instr.i[FX_DRIVE] / 32.0f;
            auto panAmt = instr.i[FX_PAN_AMT] / 512.0f;
            auto panFreq = 6.283184f * std::pow(2, (instr.i[FX_PAN_FREQ] - 9)) / rowLen;
            auto dlyAmt = instr.i[FX_DELAY_AMT] / 255.0f;
            auto dly = instr.i[FX_DELAY_TIME] * rowLen & ~1;  // Must be an even number

            // Calculate start sample number for this row in the pattern
            rowStartSample = (p * patternLen + row) * rowLen;

            // Generate notes for this pattern row
            for (col = 0; col < 4; ++col) {
                n = cp ? instr.c[cp.value() - 1].n[row + col * patternLen] : 0;
                if (n) {
                    const auto& noteBuf = find_or_create(&noteCache, n.value(), [&]{return createNote(instr, n.value(), rowLen, &random_engine); });
                    for (j = 0, i = rowStartSample * 2; j < cint(noteBuf.size()); j++, i += 2) {
                        chnBuf[i] += noteBuf[j];
                    }
                }
            }

            // Perform effects for this pattern row
            for (j = 0; j < rowLen; j++) {
                // Dry mono-sample
                k = (rowStartSample + j) * 2;
                rsample = static_cast<float>(chnBuf[k]);

                // We only do effects if we have some sound input
                if (chnBuf[k] != 0 || filterActive) {
                    // State variable filter
                    f = fxFreq;
                    if (fxLFO) {
                        f *= oscLFO(lfoFreq * k) * lfoAmt + 0.5f;
                    }
                    f = 1.5f * std::sin(f);
                    low += f * band;
                    high = q * (rsample - band) - low;
                    band += f * high;
                    rsample = fxFilter == 3 ? band : fxFilter == 1 ? high : low;

                    // Distortion
                    if (has_dist) {
                        rsample *= dist;
                        rsample = rsample < 1 ? rsample > -1 ? osc_sin(rsample*0.25f) : -1 : 1;
                        rsample /= dist;
                    }

                    // Drive
                    rsample *= drive;

                    // Is the filter active (i.e. still audiable)?
                    filterActive = (rsample * rsample) > 1e-5f;

                    // Panning
                    t = static_cast<float>(std::sin(panFreq * k)) * panAmt + 0.5f;
                    lsample = rsample * (1 - t);
                    rsample *= t;
                } else {
                    lsample = 0;
                }

                // Delay is always done, since it does not need sound input
                if (k >= dly) {
                    // Left channel = left + right[-p] * t
                    lsample += chnBuf[k-dly+1] * dlyAmt;

                    // Right channel = right + left[-p] * t
                    rsample += chnBuf[k-dly] * dlyAmt;
                }

                // Store in stereo channel buffer (needed for the delay effect)
                chnBuf[k] = truncate(lsample);
                chnBuf[k+1] = truncate(rsample);

                // ...and add to stereo mix buffer
                mMixBuf[k] += truncate(lsample);
                mMixBuf[k+1] += truncate(rsample);
            }
        }
    }

    // Next iteration. Return progress (1.0 == done!).
    mCurrentCol++;
    return mCurrentCol / static_cast<float>(mSong.numChannels);
};


// Create a AudioBuffer from the generated audio data
#if 0
Player::createAudioBuffer(context) {
    auto buffer = context.createBuffer(2, mNumWords / 2, 44100);
    for (auto i = 0; i < 2; i ++) {
        auto data = buffer.getChannelData(i);
        for (auto j = i; j < mNumWords; j += 2) {
            data[j >> 1] = mMixBuf[j] / 65536;
        }
    }
    return buffer;
};
#endif


// Create a WAVE formatted Uint8Array from the generated audio data
std::vector<u8> Player::createWave()
{
    // Create WAVE header
    auto headerLen = 44;
    int l1 = headerLen + mNumWords * 2 - 8;
    int l2 = l1 - 36;
    auto narrow = [](int n) { return static_cast<u8>(n); };
    auto wave = std::vector<u8>
    {
        82,73,70,70,
        narrow(l1 & 255),narrow((l1 >> 8) & 255),narrow((l1 >> 16) & 255),narrow((l1 >> 24) & 255),
        87,65,86,69,102,109,116,32,16,0,0,0,1,0,2,0,
        68,172,0,0,16,177,2,0,4,0,16,0,100,97,116,97,
        narrow(l2 & 255),narrow((l2 >> 8) & 255),narrow((l2 >> 16) & 255),narrow((l2 >> 24) & 255)
    };
    wave.resize(headerLen + mNumWords * 2);

    // Append actual wave data
    for (auto i = 0, idx = headerLen; i < mNumWords; ++i) {
        // Note: We clamp here
        auto y = mMixBuf[i];
        y = y < -32767 ? -32767 : (y > 32767 ? 32767 : y);
        wave[idx++] = y & 255;
        wave[idx++] = (y >> 8) & 255;
    }

    // Return the WAVE formatted typed array
    return wave;
};


// Get n samples of wave data at time t [s]. Wave data in range [-2,2].
std::vector<float> Player::getData(float t, int n)
{
    auto i = cint(2 * std::floor(t * 44100));
    auto d = std::vector<float>(n);
    for (auto j = 0; j < 2*n; j += 1) {
        auto k = i + j;
        d[j] = t > 0 && k < cint(mMixBuf.size()) ? mMixBuf[k] / 32768.0f : 0.0f;
    }
    return d;
}

