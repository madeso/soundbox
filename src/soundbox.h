#pragma once


#include <vector>
#include <random>
#include <optional>

#include "types.h"


struct Column
{
    std::vector<std::optional<int>> n;
    std::vector<std::optional<int>> f;
};

struct Instrument
{
    std::vector<int> i;

    // todo(Gustav): for replacing later
    // int osc1_waveform;
    // int osc1_vol;
    // int osc1_semi;
    // int osc1_xenv;
    // int osc2_waveform;
    // int osc2_vol;
    // int osc2_semi;
    // int osc2_detune;
    // int osc2_xenv;
    // int noise_vol;
    // int env_attack;
    // int env_sustain;
    // int env_release;
    // int env_exp_decay;
    // int arp_chord;
    // int arp_speed;
    // int lfo_waveform;
    // int lfo_amt;
    // int lfo_freq;
    // int lfo_fx_freq;
    // int fx_filter;
    // int fx_freq;
    // int fx_resonance;
    // int fx_dist;
    // int fx_drive;
    // int fx_pan_amt;
    // int fx_pan_freq;
    // int fx_delay_amt;
    // int fx_delay_time;

    std::vector<std::optional<int>> p;    // Patterns
    std::vector<Column> c; // Columns
};


struct Song
{
    std::vector<Instrument> songData;

    int rowLen;      // In sample lengths
    int patternLen;  // Rows per pattern
    int endPattern;  // End pattern
    int numChannels; // Number of channels
};


struct Player
{
    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    Player(const Song& song);


    //--------------------------------------------------------------------------
    // Public methods
    //--------------------------------------------------------------------------

    // Generate audio data for a single track
    float generate();

    // Create a WAVE formatted Uint8Array from the generated audio data
    std::vector<u8> createWave();

    // Get n samples of wave data at time t [s]. Wave data in range [-2,2].
    std::vector<float> getData(float t, int n);

    using RandomEngine = std::mt19937;

private:
    RandomEngine random_engine;

    Song mSong;
    int mLastRow;
    int mCurrentCol;
    int mNumWords;
    std::vector<i32> mMixBuf;
};

