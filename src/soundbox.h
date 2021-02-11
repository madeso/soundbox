#pragma once


#include <vector>
#include <array>
#include <random>
#include <optional>

#include "types.h"


struct Column
{
    std::vector<std::optional<int>> n;
    std::vector<std::optional<int>> f;
};

enum
{
    OSC1_WAVEFORM,
    OSC1_VOL,
    OSC1_SEMI,
    OSC1_XENV,
    OSC2_WAVEFORM,
    OSC2_VOL,
    OSC2_SEMI,
    OSC2_DETUNE,
    OSC2_XENV,
    NOISE_VOL,
    ENV_ATTACK,
    ENV_SUSTAIN,
    ENV_RELEASE,
    ENV_EXP_DECAY,
    ARP_CHORD,
    ARP_SPEED,
    LFO_WAVEFORM,
    LFO_AMT,
    LFO_FREQ,
    LFO_FX_FREQ,
    FX_FILTER,
    FX_FREQ,
    FX_RESONANCE,
    FX_DIST,
    FX_DRIVE,
    FX_PAN_AMT,
    FX_PAN_FREQ,
    FX_DELAY_AMT,
    FX_DELAY_TIME,
    INSTUMENT_DATA_COUNT
};

struct Instrument
{
    std::array<int, INSTUMENT_DATA_COUNT> i;

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

