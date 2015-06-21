/*----------------------------------------------------------------------+
 |  filename:   ExactFFT_TEST.c                                         |
 |----------------------------------------------------------------------|
 |  version:    7.10                                                    |
 |  revision:   07/09/2014  11:41                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Тест комплексного FFT                                   |
 |----------------------------------------------------------------------*/

#include "ExactFFT.c"

/// <summary>
/// Получение указателя на файловый поток
/// <summary>
/// <param name="dirName"> Имя директории. </param>
/// <param name="fileName"> Имя файла. </param>
extern FILE * GetStreamPointer(char *dirName, char *fileName, bool toWrite);

int main(int argc, char* argv[])
{
    int i, frameWidth, polyDiv2, N, N2, depth, cosTW;
    double beta, ACH_Difference, sampFreq, trueFreq, exactFreq,
           exactFreqDiff;
    short *FFT_S_short;
    double *FFT_S, *FFT_T, *MagC, *MagL, *MagR, *ACH, *ArgC, *ArgL, *ArgR, *PhaseLR;
    int FFT_S_Offset;
    
    bool useTaperWindow, recoverAfterTaperWindow,
         useNorm, direction, usePolyphase, simpleMode, isMirror, isComplex;

    CFFT_Object *fftObj;
    CFFT_SelfTestResult selfTestResult;
    FILE *testSignalFile;

    // ***************************************************

    printf("ExactFFT TEST \"C\" 7.10, (c) TESLA, 2014");

    // ***************************************************
    // * КОНСТРУКТОР
    // ***************************************************
    frameWidth = 4096;
    beta       = 28; // 28
    polyDiv2   = 1;
    cosTW  = BLACKMAN_HARRIS_92dbPS;
    fftObj = CFFT_Constructor_Cosine(frameWidth, cosTW, polyDiv2);
    //fftObj = CFFT_Constructor_Kaiser(frameWidth, beta, polyDiv2);

    // ***************************************************
    // * САМОДИАГНОСТИКА
    // ***************************************************
    ACH_Difference = 1000;
    selfTestResult = SelfTest_RND(ACH_Difference, fftObj);
    FFT_S          = (double *)calloc((frameWidth << 1), sizeof(double));
    FFT_S_short    = (short  *)calloc((frameWidth << 1), sizeof(short));
    FFT_T          = (double *)calloc((frameWidth << 1), sizeof(double));

    // (Количество точек FFT / 2) - количество гармоник вместе с нулевой
    N  = fftObj->N;
    N2 = N >> 1;

    // Массивы результатов Фурье-анализа
    MagC    = (double *)calloc(N,  sizeof(double));
    MagL    = (double *)calloc(N2, sizeof(double));
    MagR    = (double *)calloc(N2, sizeof(double));
    ACH     = (double *)calloc(N2, sizeof(double));
    ArgC    = (double *)calloc(N,  sizeof(double));
    ArgL    = (double *)calloc(N2, sizeof(double));
    ArgR    = (double *)calloc(N2, sizeof(double));
    PhaseLR = (double *)calloc(N2, sizeof(double));

    testSignalFile = GetStreamPointer(NULL,
                     "3600_Hz_STEREO_36000_SampleRate_36_deg_65536.raw",
                     FALSE);

    if(testSignalFile == NULL)
    {
        printf("\nCan't open 3600_Hz_STEREO_36000_SampleRate_36_deg_65536.raw!");
        return 1;
    }

    fread(FFT_S_short, sizeof(short), (frameWidth << 1), testSignalFile);
    for(i = 0; i < (frameWidth << 1); ++i)
    {
        FFT_S[i] = (double)FFT_S_short[i];
    }
    DumpDouble(FFT_S, (frameWidth << 1), DUMP_NAME, "FFT_S.double");

    // Прямой прогон FFT
    useTaperWindow = TRUE;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = FALSE;
    useNorm      = TRUE;
    direction    = TRUE;
    usePolyphase = FALSE;
    isMirror     = TRUE;

    CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                 recoverAfterTaperWindow, useNorm, direction,
                 usePolyphase, fftObj);
    DumpDouble(FFT_T, (frameWidth << 1), DUMP_NAME, "FFT_T.double");

    CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                 usePolyphase, fftObj);
    DumpDouble(MagL,    N2, DUMP_NAME, "MagL.double");
    DumpDouble(MagR,    N2, DUMP_NAME, "MagR.double");
    DumpDouble(PhaseLR, N2, DUMP_NAME, "PhaseLR.double");

    CFFT_ComplexExplore(FFT_T, MagC, ArgC, usePolyphase, isMirror, fftObj);
    DumpDouble(MagC,    N,  DUMP_NAME, "MagC.double");
    DumpDouble(ArgC,    N,  DUMP_NAME, "ArgC.double");

    fclose(testSignalFile);

    // Вычисление точной частоты
    sampFreq = 36000;
    trueFreq = 3600;
    depth = 20;
    isComplex = FALSE;
    exactFreq = ExactFreqAuto(MagL, depth, sampFreq, isComplex, fftObj);
    exactFreqDiff = fabs(exactFreq - trueFreq);

    DumpDouble(&exactFreq,     1, DUMP_NAME, "exactFreq.double");
    DumpDouble(&trueFreq,      1, DUMP_NAME, "trueFreq.double");
    DumpDouble(&exactFreqDiff, 1, DUMP_NAME, "exactFreqDiff.double");

    // ***************************************************
    // * ДЕСТРУКТОР
    // ***************************************************
    SAFE_DELETE(FFT_S);
    SAFE_DELETE(FFT_S_short);
    SAFE_DELETE(FFT_T);
    SAFE_DELETE(MagC);
    SAFE_DELETE(MagL);
    SAFE_DELETE(MagR);
    SAFE_DELETE(ACH);
    SAFE_DELETE(ArgC);
    SAFE_DELETE(ArgL);
    SAFE_DELETE(ArgR);
    SAFE_DELETE(PhaseLR);

    CFFT_Destructor(fftObj);

    return 0;
}