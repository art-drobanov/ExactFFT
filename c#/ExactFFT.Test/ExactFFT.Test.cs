/*----------------------------------------------------------------------+
 |  filename:   ExactFFT_TEST.cs                                        |
 |----------------------------------------------------------------------|
 |  version:    8.40                                                    |
 |  revision:   17.11.2016  16:32                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Тест комплексного FFT                                   |
 |----------------------------------------------------------------------*/

using System;
using System.IO;
using DrAF.Utilities;

namespace DrAF.DSP
{
    class ExactFFT_test
    {
        private static void Main(string[] args)
        {
            int i, j, frameWidth, polyDiv2, windowStep, N, N2, depth;
            // double kaiserBeta;
            double ACH_Difference, sampFreq, trueFreq, exactFreq, exactFreqDiff;
            double[] FFT_S, FFT_T, MagC, MagL, MagR, ACH, ArgC, ArgL, ArgR, PhaseLR;
            int FFT_S_Offset;

            bool useTaperWindow, recoverAfterTaperWindow, useNorm, direction, usePolyphase, isMirror, isComplex;

            ExactFFT.CFFT_Object fftObj;
            ExactFFT.CFFT_SelfTestResult selfTestResult;
            ExactFFT.IsDumpMode = true;

            // ***************************************************

            Console.WriteLine("ExactFFT TEST \"C#\" 8.40, (c) DrAF, 2016");

            // ***************************************************
            // * КОНСТРУКТОР
            // ***************************************************
            frameWidth = 4096;
            //kaiserBeta = 28; // 28
            polyDiv2 = 1;
            windowStep = frameWidth / 3;
            isComplex = false;
            ExactFFT.TaperWindow taperWindow = ExactFFT.TaperWindow.BLACKMAN_HARRIS_92dbPS;
            fftObj = ExactFFT.CFFT_Constructor_Cosine(frameWidth, taperWindow, polyDiv2, windowStep, isComplex);
            //fftObj = ExactFFT.CFFT_Constructor_Kaiser(frameWidth, beta, polyDiv2, windowStep, isComplex);

            // ***************************************************
            // * САМОДИАГНОСТИКА
            // ***************************************************
            ACH_Difference = 1000;
            selfTestResult = ExactFFT.SelfTest_RND(ACH_Difference, fftObj);
            Console.WriteLine("Process & Explore: {0} ms / {1} ms", (selfTestResult.CFFT_Process_time * 1000).ToString("F4"), (selfTestResult.CFFT_Explore_time * 1000).ToString("F4"));
            Console.WriteLine("Self-test result: {0}", selfTestResult.AllOK);

            // Источник и приемник
            FFT_S = new double[frameWidth << 1];
            FFT_T = new double[frameWidth << 1];

            // (Количество точек FFT / 2) - количество гармоник вместе с нулевой
            N = fftObj.N;
            N2 = N >> 1;

            // Массивы результатов Фурье-анализа
            MagC = new double[N];
            MagL = new double[N2];
            MagR = new double[N2];
            ACH  = new double[N2];
            ArgC = new double[N];
            ArgL = new double[N2];
            ArgR = new double[N2];
            PhaseLR = new double[N2];

            // Читаем все данные из файла с тестовым сигналом
            if(!File.Exists("3600_Hz_STEREO_36000_SampleRate_36_deg_65536.raw"))
            {
                Console.WriteLine("\nCan't open 3600_Hz_STEREO_36000_SampleRate_36_deg_65536.raw!");
                return;
            }

            byte[] testSignalBytes = File.ReadAllBytes("3600_Hz_STEREO_36000_SampleRate_36_deg_65536.raw");
            for(i = 0, j = 0; i < (frameWidth << 1); ++i, j += 2)
            {
                FFT_S[i] = (double)BitConverter.ToInt16(testSignalBytes, j);
            }
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "FFT_S.double", FFT_S);

            // Прямой прогон FFT
            useTaperWindow = true;
            FFT_S_Offset = 0;
            recoverAfterTaperWindow = false;
            useNorm      = true;
            direction    = true;
            usePolyphase = false;
            isMirror     = true;

            ExactFFT.CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                                  recoverAfterTaperWindow, useNorm, direction,
                                  usePolyphase, fftObj);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "FFT_T.double", FFT_T);

            ExactFFT.CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                                  usePolyphase, fftObj);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "MagL.double", MagL);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "MagR.double", MagR);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "PhaseLR.double", PhaseLR);

            ExactFFT.CFFT_ComplexExplore(FFT_T, MagC, ArgC, usePolyphase, isMirror, fftObj);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "MagC.double", MagC);
            DebugHelper.WriteDoubles(ExactFFT.DumpName, "ArgC.double", ArgC);

            // Вычисление точной частоты
            sampFreq = 36000;
            trueFreq = 3600;
            depth = 20;
            exactFreq = ExactFFT.ExactFreqAuto(MagL, depth, sampFreq, fftObj);
            exactFreqDiff = Math.Abs(exactFreq - trueFreq);

            DebugHelper.WriteDouble(ExactFFT.DumpName, "exactFreq.double", exactFreq);
            DebugHelper.WriteDouble(ExactFFT.DumpName, "trueFreq.double", trueFreq);
            DebugHelper.WriteDouble(ExactFFT.DumpName, "exactFreqDiff.double", exactFreqDiff);

            // ***************************************************
            // * ДЕСТРУКТОР
            // ***************************************************
            FFT_S   = null;
            FFT_T   = null;
            MagC    = null;
            MagL    = null;
            MagR    = null;
            ACH     = null;
            ArgC    = null;
            ArgL    = null;
            ArgR    = null;
            PhaseLR = null;

            ExactFFT.CFFT_Destructor(fftObj);
        }
    }
}
