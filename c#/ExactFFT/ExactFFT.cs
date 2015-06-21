/*----------------------------------------------------------------------+
 |  filename:   ExactFFT.cs                                             |
 |----------------------------------------------------------------------|
 |  version:    8.00                                                    |
 |  revision:   25/11/2014  17:43                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Комплексное FFT                                         |
 |----------------------------------------------------------------------*/

using System;
using System.IO;
using DrAF.Utilities;

namespace DrAF.DSP
{
    /// <summary>
    /// Комплексное FFT
    /// </summary>
    public static class ExactFFT
    {
        /// <summary> Имя дампа (имя директории). </summary>
        public static string DumpName = ".NET.dump";

        /// <summary> Активирован режим дампа? </summary>
        public static bool IsDumpMode = false;

        // Константы
        //                                          3.14159265358979323846264338328
        public const double M_PI            = 3.14159265358979323846; //..264338328
        public const double M_2PI           = 2 * M_PI;
        public const double FLOAT_MIN       = 3.4E-38;    // Допустимый минимум для операций float
        public const double MAX_FFT_DIFF    = 1E-7;       // Максимальная погрешность FFT
        public const double MIN_FRAME_WIDTH = 8;          // Наименьший "рабочий" размер окна FFT
        public const double MAX_KAISER_BETA = 28;         // Max. Beta (Kaiser window), SL: ~240 dB
        public const double MAX_PATH        = 256;        // Максимальная длина пути        

        // "Булевские" константы
        public const bool DIRECT                 = true;  // Обозначение прямого прохода FFT
        public const bool REVERSE                = false; // Обозначение обратного прохода FFT
        public const bool USING_NORM             = true;  // Используется нормализация
        public const bool NOT_USING_NORM         = false; // Не используется нормализация
        public const bool USING_TAPER_WINDOW     = true;  // Используется взвешивающее окно
        public const bool NOT_USING_TAPER_WINDOW = false; // Не используется взвешивающее окно
        public const bool USING_POLYPHASE        = true;  // Используется полифазное FFT
        public const bool NOT_USING_POLYPHASE    = false; // Не используется полифазное FFT

        /// <summary>
        /// Структура "Объект FFT"
        /// </summary>
        public class CFFT_Object
        {
            //-------------------------------------------------------------------------
            public int N;           // Количество точек FFT
            public int NN;          // Кол-во чисел (re + im) FFT
            public int NPoly;       // Пересчитанное для полифазного FFT количество точек
            public int NNPoly;      // Кол-во чисел (re + im) полифазного FFT
            public double Beta;     // Формирующий коэффициент "Beta" окна Кайзера
            public CosTW CosTW;     // Тип косинусного взвешивающего окна (если нет - используется окно Кайзера)
            public int PolyDiv;     // Делитель "полифазности" FFT ("0" - обычное FFT)
            public int WindowStep;  // Шаг окна FFT
            public bool IsComplex;  // Используется комплексный вход? (бывает COMPLEX и L+R)
            //-------------------------------------------------------------------------
            public int[] FFT_P;     // Вектор изменения порядка следования данных перед FFT
            public int[] FFT_PP;    // Вектор изменения порядка... (для полифазного FFT)
            public double[] FFT_TW; // Взвешивающее окно
            //-------------------------------------------------------------------------
        }

        /// <summary>
        /// Результат самодиагностики
        /// </summary>
        public struct CFFT_SelfTestResult
        {
            //-------------------------------------------------------------------------
            public int AllOK; // Результат тестирования точности
            //-------------------------------------------------------------------------
            public double MaxDiff_ACH; // Максимальная невязка по расчету заданной АЧХ
            public double MaxDiff_ALG_to_EXP_to_ALG; // Max. невязка ALG . EXP и обратно
            public double MaxDiff_FORWARD_BACKWARD;  // Max. невязка FORVARD + BACKWARD
            public double MaxDiff_FORWARD_BACKWARD_AntiTW; //...то же + восст. после TW
            public double MaxDiff_PhaseLR;   // Макс. невязка по расчету разности хода фаз
            public double CFFT_Process_time; // Время работы CFFT_Process()
            public double CFFT_Explore_time; // Время работы CFFT_Explore()
            //-------------------------------------------------------------------------
        }
      
        /// <summary>
        /// Логарифм по произвольному основанию
        /// </summary>
        /// <param name="arg"> Аргумент логарифма. </param>
        /// <param name="logBase"> Основание логарифма. </param>
        public static double LogX(double arg, double logBase)
        {
            return Math.Log(arg) / Math.Log(logBase);
        }

        /// <summary>
        /// Приведение значения к ближайшей снизу степени двойки
        /// </summary>
        /// <param name="arg"> Входной аргумент. </param>
        public static int ToLowerPowerOf2(int arg)
        {
            return (int)Math.Pow(2, (int)(LogX(arg, 2)));
        }

        /// <summary>
        /// Переход от линейной шкалы к dB
        /// </summary>
        /// <param name="arg"> Входной аргумент (линейная шкала). </param>
        /// <param name="zero_db_level"> "Нулевой" уровень dB. </param>
        public static double To_dB(double arg, double zero_db_level)
        {
            return 10.0 * Math.Log(arg / zero_db_level); // log
        }

        /// <summary>
        /// Переход от dB к линейной шкале
        /// </summary>
        /// <param name="arg"> Входной аргумент (логарифмическая шкала). </param>
        /// <param name="zero_db_level"> "Нулевой" уровень dB. </param>
        public static double From_dB(double arg, double zero_db_level)
        {
            return zero_db_level * Math.Pow(10, arg / 10.0); // exp
        }

        /// <summary>
        /// Получение частоты заданного узла FFT
        /// </summary>
        /// <param name="FFT_Node"> Номер гармоники. </param>
        /// <param name="sampFreq"> Частота семплирования. </param>        
        /// <param name="N"> Размер окна FFT. </param>
        /// <param name="isComplex"> Комплексный режим? </param>
        public static double FreqNode(double FFT_Node, double sampFreq, int N, bool isComplex)
        {
            return (FFT_Node * sampFreq) / ((double)N * (isComplex ? 2.0 : 1.0));
        }

        /// <summary>
        /// Получение узла FFT по заданной частоте
        /// </summary>
        /// <param name="freqNode"> Заданная частота. </param>
        /// <param name="sampFreq"> Частота семплирования. </param>        
        /// <param name="N"> Размер окна FFT. </param>
        /// <param name="isComplex"> Комплексный режим? </param>
        public static double FFT_Node(double freqNode, double sampFreq, int N, bool isComplex)
        {
            return (freqNode * ((double)N * (isComplex ? 2.0 : 1.0))) / sampFreq;
        }

        /// <summary>
        /// Нормирование фазы
        /// </summary>
        /// <param name="phase"> Значение фазы для нормирования. </param>
        public static double PhaseNorm(double phase)
        {
            if(phase > 0) while(phase >= M_PI)  phase -= M_2PI;
                     else while(phase <= -M_PI) phase += M_2PI;

            return phase;
        }

        /// <summary>
        /// Безопасное в смысле значений аргументов вычисление арктангенса
        /// </summary>
        /// <param name="im"> Мнимая часть комплексного числа. </param>
        /// <param name="re"> Действительная часть комплексного числа. </param>
        public static double Safe_atan2(double im, double re)
        {
            return (((re < 0) ? -re : re) < FLOAT_MIN) ? 0 : Math.Atan2(im, re);
        }

        /// <summary>
        /// Заполнение вектора изменения порядка следования данных перед FFT
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void fill_FFT_P(CFFT_Object fftObj)
        {
            int i, j, shift;

            // Выделяем память под вектор перестановок FFT...            
            fftObj.FFT_P = new int[fftObj.NN];

            // Заполняем вектор изменения порядка следования данных...
            for(j = 0; j < LogX(fftObj.N, 2); ++j)
            {
                for(i = 0; i < fftObj.N; ++i)
                {
                    fftObj.FFT_P[i << 1] = ((fftObj.FFT_P[i << 1] << 1) +
                                            ((i >> j) & 1));
                }
            }

            shift = (fftObj.FFT_P[2] == (fftObj.N >> 1)) ? 1 : 0;
            for(i = 0; i < fftObj.NN; i += 2)
            {
                fftObj.FFT_P[i + 1] = (fftObj.FFT_P[i + 0] <<= shift) + 1;
            }
        }

        /// <summary>
        /// Заполнение вектора изменения порядка следования данных перед FFT
        /// (для полифазного FFT)
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void fill_FFT_PP(CFFT_Object fftObj)
        {
            int i, j;

            // Выделяем память под вектор перестановок FFT...
            fftObj.FFT_PP = new int[fftObj.NNPoly];

            // Заполняем вектор изменения порядка следования данных
            // (для полифазного FFT)...
            for(j = 0; j < LogX(fftObj.NPoly, 2); ++j)
            {
                for(i = 0; i < fftObj.NPoly; ++i)
                {
                    fftObj.FFT_PP[i << 1] = ((fftObj.FFT_PP[i << 1] << 1) +
                                             ((i >> j) & 1));
                }
            }

            for(i = 0; i < fftObj.NNPoly; i += 2)
            {
                fftObj.FFT_PP[i + 1] = (fftObj.FFT_PP[i + 0] <<= 1) + 1;
            }
        }

        /// <summary>
        /// Типы взвешивающих окон FFT
        /// </summary>
        /// Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)        
        public enum CosTW
        {
            NONE,
            RECTANGULAR_13dbPS,
            HANN_31dbPS,
            HAMMING_43dbPS,
            MAX_ROLLOFF_3_TERM_46dbPS,
            BLACKMAN_58dbPS,
            COMPROMISE_3_TERM_64dbPS,
            EXACT_BLACKMAN_68dbPS,
            MIN_SIDELOBE_3_TERM_71dbPS,
            MAX_ROLLOFF_4_TERM_60dbPS,
            COMPROMISE1_4_TERM_82dbPS,
            COMPROMISE2_4_TERM_93dbPS,
            BLACKMAN_HARRIS_92dbPS,
            NUTTALL_93dbPS,
            BLACKMAN_NUTTALL_98dbPS,
            ROSENFIELD
        };

        /// <summary>
        /// Заполнение вектора косинусного взвешивающего окна
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void fill_FFT_TW_Cosine(CFFT_Object fftObj)
        {
            // Выделяем память под взвешивающее окно...
            fftObj.FFT_TW = new double[fftObj.NN];

            // Формирующие параметры взвешивающего косинусного окна
            double a0, a1, a2, a3, ad;

            // Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)
            switch(fftObj.CosTW)
            {                
                case CosTW.RECTANGULAR_13dbPS:         { a0 = 1.0;       a1 = 0;         a2 = 0;         a3 = 0;         ad = 1.0;     break; }
                case CosTW.HANN_31dbPS:                { a0 = 1.0;       a1 = 1.0;       a2 = 0;         a3 = 0;         ad = 2;       break; }
                case CosTW.HAMMING_43dbPS:             { a0 = 0.54;      a1 = 0.46;      a2 = 0;         a3 = 0;         ad = 1.0;     break; }
                case CosTW.MAX_ROLLOFF_3_TERM_46dbPS:  { a0 = 0.375;     a1 = 0.5;       a2 = 0.125;     a3 = 0;         ad = 1.0;     break; }
                case CosTW.BLACKMAN_58dbPS:            { a0 = 0.42;      a1 = 0.5;       a2 = 0.08;      a3 = 0;         ad = 1.0;     break; }
                case CosTW.COMPROMISE_3_TERM_64dbPS:   { a0 = 0.40897;   a1 = 0.5;       a2 = 0.09103;   a3 = 0;         ad = 1.0;     break; }
                case CosTW.EXACT_BLACKMAN_68dbPS:      { a0 = 7938.0;    a1 = 9240.0;    a2 = 1430.0;    a3 = 0;         ad = 18608.0; break; }
                case CosTW.MIN_SIDELOBE_3_TERM_71dbPS: { a0 = 0.4243801; a1 = 0.4973406; a2 = 0.0782793; a3 = 0;         ad = 1.0;     break; }
                case CosTW.MAX_ROLLOFF_4_TERM_60dbPS:  { a0 = 10.0;      a1 = 15.0;      a2 = 6.0;       a3 = 1;         ad = 32.0;    break; }
                case CosTW.COMPROMISE1_4_TERM_82dbPS:  { a0 = 0.338946;  a1 = 0.481973;  a2 = 0.161054;  a3 = 0.018027;  ad = 1.0;     break; }
                case CosTW.COMPROMISE2_4_TERM_93dbPS:  { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
                default:
                case CosTW.BLACKMAN_HARRIS_92dbPS:     { a0 = 0.35875;   a1 = 0.48829;   a2 = 0.14128;   a3 = 0.01168;   ad = 1.0;     break; }
                case CosTW.NUTTALL_93dbPS:             { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
                case CosTW.BLACKMAN_NUTTALL_98dbPS:    { a0 = 0.3635819; a1 = 0.4891775; a2 = 0.1365995; a3 = 0.0106411; ad = 1.0;     break; }
                case CosTW.ROSENFIELD:                 { a0 = 0.762;     a1 = 1.0;       a2 = 0.238;     a3 = 0;         ad = a0;      break; }
            }

            // Заполняем взвешивающее окно коэффициентами...
            for(int i = 0; i < fftObj.N; ++i)
            {
                double arg  = (2.0 * Math.PI * i) / (double)fftObj.N;
                double wval = (a0 - a1 * Math.Cos(arg) + a2 * Math.Cos(2 * arg) - a3 * Math.Cos(3 * arg)) / ad;
                fftObj.FFT_TW[(i << 1) + 1] = fftObj.FFT_TW[(i << 1) + 0] = wval;
            }
        }

        /// <summary>
        /// Модифицир. функция Бесселя нулевого порядка первого рода
        /// </summary>
        /// <param name="arg"> Аргумент функции. </param>
        /// <returns> Значение функции. </returns>
        public static double BesselI0(double arg)
        {
            double numerator, denominator, z, z1, z2, z3, z4, z5,
                   z6, z7, z8, z9, z10, z11, z12, z13, z_1, z_2;

            if(arg == 0.0)
            {
                return 1.0;
            }

            z = arg * arg;

            z1  = z * 0.210580722890567e-22 + 0.380715242345326e-19;
            z2  = z * z1 + 0.479440257548300e-16;
            z3  = z * z2 + 0.435125971262668e-13;
            z4  = z * z3 + 0.300931127112960e-10;
            z5  = z * z4 + 0.160224679395361e-7;
            z6  = z * z5 + 0.654858370096785e-5;
            z7  = z * z6 + 0.202591084143397e-2;
            z8  = z * z7 + 0.463076284721000e0;
            z9  = z * z8 + 0.754337328948189e2;
            z10 = z * z9 + 0.830792541809429e4;
            z11 = z * z10 + 0.571661130563785e6;
            z12 = z * z11 + 0.216415572361227e8;
            z13 = z * z12 + 0.356644482244025e9;

            numerator = z * z13 + 0.144048298227235e10;

            z_1 = z - 0.307646912682801e4;
            z_2 = z * z_1 + 0.347626332405882e7;
            denominator = z * z_2 - 0.144048298227235e10;

            return -numerator / denominator;
        }

        /// <summary>
        /// Метод заполнения взвешивающего окна (окно Кайзера)
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void fill_FFT_TW_Kaiser(CFFT_Object fftObj)
        {
            int i, j;
            double norm, arg, w;

            // Выделяем память под вектор перестановок FFT...
            fftObj.FFT_TW = new double[fftObj.NN];

            // Нормирующий коэффициент окна Кайзера
            norm = BesselI0(fftObj.Beta);

            // Заполняем взвешивающее окно...
            for(i = 1; i <= (fftObj.N >> 1); ++i)
            {
                // arg = Beta * sqrt(1-(((2*(i-1))/(N-1))-1)^2);
                arg = fftObj.Beta *
                                    Math.Sqrt(
                                               1 - Math.Pow(
                                                             (
                                                                (double)((i - 1) << 1)
                                                              /
                                                                (double)(fftObj.N - 1)
                                                             ) - 1
                                                        , 2)
                                              );
                
                w = BesselI0(arg) / norm;

                j = i - 1; // Приводим индекс от базы "1" к базе "0"
                fftObj.FFT_TW[(j << 1) + 0] = w; // left re
                fftObj.FFT_TW[(j << 1) + 1] = w; // left im
                fftObj.FFT_TW[(fftObj.NN - 2) - (j << 1) + 0] = w; // right re
                fftObj.FFT_TW[(fftObj.NN - 2) - (j << 1) + 1] = w; // right im
            }
        }

        /// <summary>
        /// "Контролер" объекта FFT
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Булевский флаг корректности конфигурации FFT. </returns>
        public static bool CFFT_Inspector(CFFT_Object fftObj)
        {
            // Размер окна не может быть меньше предельно-допустимого!
            if((fftObj.N     < MIN_FRAME_WIDTH) ||
               (fftObj.NPoly < MIN_FRAME_WIDTH) ||
               (fftObj.Beta  > MAX_KAISER_BETA))
            {
                return false;
            }

            return true;
        }

        /// <summary>
        /// "Деструктор" объекта FFT
        /// </summary>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void CFFT_Destructor(CFFT_Object fftObj)
        {
            fftObj.FFT_P  = null;
            fftObj.FFT_PP = null;
            fftObj.FFT_TW = null;
        }

        /// <summary>
        /// Фабрика объектов FFT
        /// </summary>
        /// <param name="frameWidth"> Размер кадра. </param>
        /// <param name="cosTW"> Тип косинусного взвешивающего окна (если нет - используется окно Кайзера). </param>
        /// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
        /// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
        /// <param name="windowStep"> Шаг окна FFT. </param>
        /// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
        public static CFFT_Object CFFT_Init(int frameWidth, CosTW cosTW, double beta, int polyDiv2,
                                            int windowStep, bool isComplex)
        {
            // Объект-результат
            CFFT_Object fftObj = new CFFT_Object();

            // Заполнение полей объекта
            fftObj.N = ToLowerPowerOf2(frameWidth); // Размер кадра FFT
            fftObj.NN = fftObj.N << 1;              // Кол-во точек (re + im)
            fftObj.NPoly = fftObj.N >> polyDiv2;    // Размер полифазного кадра FFT
            fftObj.NNPoly = fftObj.NPoly << 1;      // Кол-во точек полифазного FFT
            fftObj.CosTW = cosTW;                   // Тип косинусного взвешивающего окна
            fftObj.Beta = beta;                     // Форм-ий коэфф. окна Кайзера
            fftObj.PolyDiv = 1 << polyDiv2;         // Полифазный делитель
            fftObj.WindowStep = windowStep;         // Шаг окна FFT
            fftObj.IsComplex = isComplex;           // Используется комплексный вход? (бывает COMPLEX и L+R)

            fill_FFT_P(fftObj);  // Вектор изменения порядка след. данных перед FFT
            fill_FFT_PP(fftObj); // Вектор изменения порядка... (для полифазного FFT)

            if(fftObj.CosTW == CosTW.NONE) //...если не задано взвешивающее окно косинусного типа
            {
                fill_FFT_TW_Kaiser(fftObj); // Взвешивающее окно Кайзера

            } else
            {
                fill_FFT_TW_Cosine(fftObj); // Косинусное взвешивающее окно
            }

            // Обрабатываем ситуацию со сбросом дампа...
            if(IsDumpMode)
            {
                Directory.CreateDirectory(DumpName);

                DebugHelper.WriteInts(DumpName,    "FFT_P.int32",   fftObj.FFT_P);
                DebugHelper.WriteInts(DumpName,    "FFT_PP.int32",  fftObj.FFT_PP);
                DebugHelper.WriteDoubles(DumpName, "FFT_TW.double", fftObj.FFT_TW);
            }

            // Если некоторые параметры не соответствуют норме...
            if(!CFFT_Inspector(fftObj))
            {
                //...- убираем объект.
                CFFT_Destructor(fftObj);
            }

            // Возвращаем объект "FFT"
            return fftObj;
        }

        /// <summary>
        /// Фабрика объектов FFT
        /// </summary>
        /// <param name="frameWidth"> Размер кадра. </param>
        /// <param name="cosTW"> Тип косинусного взвешивающего окна. </param>
        /// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
        /// <param name="windowStep"> Шаг окна FFT. </param>
        /// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
        public static CFFT_Object CFFT_Constructor_Cosine(int frameWidth, CosTW cosTW, int polyDiv2,
                                                          int windowStep, bool isComplex)
        {
            // Возвращаем объект "FFT"
            return CFFT_Init(frameWidth, cosTW, MAX_KAISER_BETA, polyDiv2, windowStep, isComplex);
        }
    
        /// <summary>
        /// Фабрика объектов FFT
        /// </summary>
        /// <param name="frameWidth"> Размер кадра. </param>
        /// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
        /// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
        /// <param name="windowStep"> Шаг окна FFT. </param>
        /// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
        public static CFFT_Object CFFT_Constructor_Kaiser(int frameWidth, double beta, int polyDiv2,
                                                          int windowStep, bool isComplex)
        {
            // Возвращаем объект "FFT"
            return CFFT_Init(frameWidth, CosTW.NONE, beta, polyDiv2, windowStep, isComplex);
        }

        /// <summary>
        /// Основной метод комплексного FFT
        /// </summary>
        /// <param name="FFT_S"> Вектор входных данных
        /// ("левый" и "правый" каналы - чет./нечет.). </param>
        /// <param name="FFT_S_Offset"> Смещение данных для анализа во
        /// входном векторе FFT_S. </param>
        /// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
        /// <param name="useTaperWindow"> Использовать взвешивающее окно? </param>
        /// <param name="recoverAfterTaperWindow"> Аннигилировать действие
        /// взвешивающего окна на обратном проходе? </param>
        /// <param name="useNorm"> Использовать нормализацию 1/N? </param>
        /// <param name="direction"> Направление преобразования (true - прямое).
        /// </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void CFFT_Process(double[] FFT_S, int FFT_S_Offset, double[] FFT_T,
                                        bool useTaperWindow, bool recoverAfterTaperWindow,
                                        bool useNorm, bool direction, bool usePolyphase,
                                        CFFT_Object fftObj)
        {
            int i, j, mmax, isteps, NN, istep, ii, m, jj;
            double isign, theta, sin05Th, wpr, wpi, wr, wi, tempr, tempi, wtemp;

            // Использование взвешивающего окна допустимо
            // только при прямом преобразовании
            if(direction && useTaperWindow)
            {
                if(!usePolyphase)
                {
                    // Обычное FFT
                    // только тогда, когда оно активно
                    for(i = 0; i < fftObj.NN; ++i)
                    {
                        FFT_T[i] = fftObj.FFT_TW[fftObj.FFT_P[i]] *
                                   FFT_S[fftObj.FFT_P[i] + FFT_S_Offset];
                    }
                }
                else
                {
                    // Полифазное FFT
                    // предполагает применение только на прямом проходе
                    for(i = 0; i < fftObj.NNPoly; ++i)
                    {
                        FFT_T[i] = 0;

                        // Накапливаем сумму текущей точки (в соответствии
                        // с количеством сегментов)
                        for(j = 0; j < fftObj.PolyDiv; ++j)
                        {
                            FFT_T[i] += fftObj.FFT_TW[fftObj.FFT_PP[i] +
                                        (j * fftObj.NNPoly)] * FFT_S[fftObj.FFT_PP[i] +
                                        (j * fftObj.NNPoly)  + FFT_S_Offset];
                        }
                    }
                }
            }
            else
            {
                // Обратный проход или прямой...
                // но без взвешивающего окна
                for(i = 0; i < fftObj.NN; ++i)
                {
                    FFT_T[i] = FFT_S[fftObj.FFT_P[i] + FFT_S_Offset];
                }
            }

            // Нормализация коэффициентов производится при её выборе и только на
            // прямом проходе алгоритма (или если ситуация 100% симметрична)
            if((!direction) && (!useNorm))
            {
                for(i = 0; i < fftObj.NNPoly; ++i)
                {
                    FFT_T[i] /= fftObj.N;
                }
            }

            // FFT Routine
            isign  = direction ? -1 : 1;
            mmax   = 2;
            isteps = 1;
            NN = usePolyphase ? fftObj.NNPoly : fftObj.NN;
            while(NN > mmax)
            {
                isteps++;
                istep   = mmax << 1;
                theta   = isign * ((2 * M_PI) / mmax);
                sin05Th = Math.Sin(0.5 * theta);
                wpr = -(2.0 * (sin05Th * sin05Th));
                wpi = Math.Sin(theta);
                wr  = 1.0;
                wi  = 0.0;

                for(ii = 1; ii <= (mmax >> 1); ++ii)
                {
                    m = (ii << 1) - 1;
                    for(jj = 0; jj <= ((NN - m) >> isteps); ++jj)
                    {
                        i = m + (jj << isteps);
                        j = i + mmax;
                        tempr = wr * FFT_T[j - 1] - wi * FFT_T[j];
                        tempi = wi * FFT_T[j - 1] + wr * FFT_T[j];
                        FFT_T[j - 1]  = FFT_T[i - 1] - tempr;
                        FFT_T[j - 0]  = FFT_T[i - 0] - tempi;
                        FFT_T[i - 1] += tempr;
                        FFT_T[i - 0] += tempi;
                    }
                    wtemp = wr;
                    wr = wr * wpr - wi * wpi + wr;
                    wi = wi * wpr + wtemp * wpi + wi;
                }
                mmax = istep;
            }

            // Нормализация коэффициентов производится при её выборе и только
            // на прямом проходе алгоритма (или если ситуация 100% симметрична)
            if(direction && useNorm)
            {
                for(i = 0; i < NN; ++i)
                {
                    FFT_T[i] /= fftObj.N;
                }
            }

            // Аннигилируем взвешивающее окно (если оно равно нулю в некоторых точках
            // - результат восстановления неизвестен и полагаем его равным нулю)
            if((!direction) && useTaperWindow && recoverAfterTaperWindow)
            {
                for(i = 0; i < fftObj.NN; ++i)
                {
                    FFT_T[i] = ((fftObj.FFT_TW[i] == 0) ?
                               0 : (FFT_T[i] / fftObj.FFT_TW[i]));
                }
            }
        }

        /// <summary>
        /// Исследование "левого" и "правого" каналов: ("левый" -
        /// действительная часть исходных данных, "правый" - мнимая часть)
        /// в режиме (L+R)
        /// </summary>
        /// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
        /// <param name="MagL"> Магнитуды "левого" канала. </param>
        /// <param name="MagR"> Магнитуды "правого" канала. </param>
        /// <param name="ACH"> АЧХ (отношение магнитуды "правого" канала к магнитуде
        /// "левого" - как "выход" / "вход"). </param>
        /// <param name="ArgL"> Аргумент "левого" канала. </param>
        /// <param name="ArgR"> Аргумент "правого" канала. </param>
        /// <param name="PhaseLR"> Разность хода фаз каналов ("правый" минус
        /// "левый"). </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void CFFT_Explore(double[] FFT_T, double[] MagL, double[] MagR,
                                        double[] ACH, double[] ArgL, double[] ArgR,
                                        double[] PhaseLR, bool usePolyphase,
                                        CFFT_Object fftObj)
        {
            int i, N;
            double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
                   lx, ly, rx, ry, argL, argR;

            // Определяем размерность FFT
            N = usePolyphase ? fftObj.NPoly : fftObj.N;

            // Вычисляем значения искомых величин для начальной точки
            // ("нулевая" гармоника)
            magL = FFT_T[0];
            magR = FFT_T[1];
            if(MagL    != null) MagL[0]    = magL;
            if(MagR    != null) MagR[0]    = magR;
            if(ACH     != null) ACH[0]     = magR / ((magL == 0) ? FLOAT_MIN : magL);
            if(ArgL    != null) ArgL[0]    = M_PI;
            if(ArgR    != null) ArgR[0]    = M_PI;
            if(PhaseLR != null) PhaseLR[0] = 0;
           
            // Работа с гармоническими точками
            for(i = 1; i < (N >> 1); ++i)
            {
                FFT_T_i_Re   = FFT_T[(i << 1) + 0];
                FFT_T_i_Im   = FFT_T[(i << 1) + 1];
                FFT_T_N_i_Re = FFT_T[((N - i) << 1) + 0];
                FFT_T_N_i_Im = FFT_T[((N - i) << 1) + 1];

                lx = FFT_T_i_Re   + FFT_T_N_i_Re;
                ly = FFT_T_i_Im   - FFT_T_N_i_Im;
                rx = FFT_T_i_Im   + FFT_T_N_i_Im;
                ry = FFT_T_N_i_Re - FFT_T_i_Re;

                magL = Math.Sqrt((lx * lx) + (ly * ly)) * 0.5;
                magR = Math.Sqrt((rx * rx) + (ry * ry)) * 0.5;
                argL = Safe_atan2(ly, lx);
                argR = Safe_atan2(ry, rx);

                if(MagL    != null) MagL[i]    = magL;
                if(MagR    != null) MagR[i]    = magR;
                if(ACH     != null) ACH[i]     = magR / ((magL == 0) ? FLOAT_MIN : magL);
                if(ArgL    != null) ArgL[i]    = argL;
                if(ArgR    != null) ArgR[i]    = argR;
                if(PhaseLR != null) PhaseLR[i] = PhaseNorm(argR - argL);
            }
        }

        /// <summary>
        /// Исследование результатов комплексного FFT (идентично CFFT из MathCAD)
        /// </summary>
        /// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
        /// <param name="Mag"> Магнитуды. </param>
        /// <param name="Arg"> Аргументы. </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="isMirror"> Зеркальное отображение спектра? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void CFFT_ComplexExplore(double[] FFT_T, double[] Mag, double[] Arg,
                                               bool usePolyphase, bool isMirror,
                                               CFFT_Object fftObj)
        {
            int i, N, N_2;
            double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
                   lx, ly, rx, ry, argL, argR;

            // Определяем размерность FFT
            N   = usePolyphase ? fftObj.NPoly : fftObj.N;
            N_2 = N >> 1;

            // Вычисляем значения искомых величин для начальной точки
            // ("нулевая" гармоника)
            lx = FFT_T[0];
            ly = FFT_T[1];
            rx = FFT_T[N + 0];
            ry = FFT_T[N + 1];
                        
            if(Mag != null)
            {
                magL = Math.Sqrt((lx * lx) + (ly * ly));
                magR = Math.Sqrt((rx * rx) + (ry * ry));
                
                if(!isMirror)
                {
                    Mag[0]   = magL;
                    Mag[N_2] = magR;

                } else
                {
                    Mag[N_2 - 1] = magL;
                    Mag[N   - 1] = magR;
                }
            }
            
            if(Arg != null)
            {
                argL = Safe_atan2(ly, lx);
                argR = Safe_atan2(ry, rx);
                
                if(!isMirror)
                {
                    Arg[0]   = argL;
                    Arg[N_2] = argR;

                } else
                {
                    Arg[N_2 - 1] = argL;
                    Arg[N   - 1] = argR;
                }
            }
            
            // Работа с гармоническими точками
            for(i = 1; i < N_2; ++i)
            {
                FFT_T_i_Re   = FFT_T[(i << 1) + 0];
                FFT_T_i_Im   = FFT_T[(i << 1) + 1];
                FFT_T_N_i_Re = FFT_T[((N - i) << 1) + 0];
                FFT_T_N_i_Im = FFT_T[((N - i) << 1) + 1];

                lx = FFT_T_i_Re;
                ly = FFT_T_i_Im;
                rx = FFT_T_N_i_Re;
                ry = FFT_T_N_i_Im;
                                                
                if(Mag != null)
                {
                    magL = Math.Sqrt((lx * lx) + (ly * ly));
                    magR = Math.Sqrt((rx * rx) + (ry * ry));

                    if(!isMirror)
                    {
                        Mag[i]     = magL;
                        Mag[N - i] = magR;

                    } else
                    {
                        Mag[N_2 - i - 1] = magL;
                        Mag[N_2 + i - 1] = magR;
                    }
                }

                if(Arg != null)
                {
                    argL = Safe_atan2(ly, lx);
                    argR = Safe_atan2(ry, rx);

                    if(!isMirror)
                    {
                        Arg[i]     = argL;
                        Arg[N - i] = argR;

                    } else
                    {
                        Arg[N_2 - i - 1] = argL;
                        Arg[N_2 + i - 1] = argR;
                    }
                }                
            }
        }

        /// <summary>
        /// Перевод значений массива double в форму dB
        /// </summary>
        /// <param name="data"> Данные для обработки. </param>
        /// <param name="zero_db_level"> Уровень "нулевого" уровня. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void dB_Scale(double[] Mag, double zero_db_level,
                                    CFFT_Object fftObj)
        {
            int i;
            for(i = 0; i < (fftObj.N * ((fftObj.IsComplex ? 2 : 1))) >> 1; ++i)
            {
                Mag[i] = 10.0 * Math.Log(Mag[i] / zero_db_level); // log
            }
        }

        /// <summary>
        /// Самотестирование внутренней точности в цикле прямого-обратного
        /// преобразования на пользовательских данных
        /// </summary>
        /// <param name="FFT_S"> Вектор входных данных ("левый" и "правый" каналы
        /// - чет./нечет.) </param>
        /// <param name="ACH_Difference"> Коэффициент превосходства правого канала
        /// по уровню в ходе проводимого теста. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Структура "Результат тестирования внутренней точности
        /// прямого-обратного комплексного преобразования Фурье". </returns>
        public static CFFT_SelfTestResult SelfTest_S(double[] FFT_S, double ACH_Difference,
                                                     CFFT_Object fftObj)
        {
            double[] FFT_S_backward, FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR;
            int N2, FFT_S_Offset, i;
            bool useTaperWindow, recoverAfterTaperWindow, useNorm, direction, usePolyphase;
            double maxDiff, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
                   lx, ly, rx, ry, lx_, ly_, rx_, ry_, currentDiff;
            long startCounter, CFFT_Process_counter, CFFT_Explore_counter, timerFrequency;
            int N_iters = 10000;

            // Cтруктура "Результат тестирования внутренней точности
            // прямого-обратного комплексного преобразования Фурье"
            CFFT_SelfTestResult selfTestResult;

            // Массив исходных данных - для заполнения на обратном ходе FFT
            FFT_S_backward = new double[fftObj.NN];

            // Целевой массив
            FFT_T = new double[fftObj.NN];

            // (Количество точек FFT / 2) - количество гармоник вместе с нулевой
            N2 = fftObj.N >> 1;

            // Массивы результатов Фурье-анализа
            MagL    = new double[N2];
            MagR    = new double[N2];
            ACH     = new double[N2];
            ArgL    = new double[N2];
            ArgR    = new double[N2];
            PhaseLR = new double[N2];

            // Не используем взвешивающее окно, но работаем
            // с нормализацией - направление прямое
            useTaperWindow = false;
            FFT_S_Offset   = 0;
            recoverAfterTaperWindow = false;
            useNorm      = true;
            direction    = true;
            usePolyphase = false;
            CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                         recoverAfterTaperWindow, useNorm, direction,
                         usePolyphase, fftObj);

            // Пользуясь результатами прямого преобразования извлекаем
            // все возможные величины
            CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                         usePolyphase, fftObj);

            // Проверяем правильность расчета магнитуд и фаз - пытаемся получить
            // исходные комплексные числа
            maxDiff = 0;
            for(i = 1; i < N2; ++i)
            {
                FFT_T_i_Re   = FFT_T[(i << 1) + 0];
                FFT_T_i_Im   = FFT_T[(i << 1) + 1];
                FFT_T_N_i_Re = FFT_T[((fftObj.N - i) << 1) + 0];
                FFT_T_N_i_Im = FFT_T[((fftObj.N - i) << 1) + 1];

                lx = FFT_T_i_Re   + FFT_T_N_i_Re;
                ly = FFT_T_i_Im   - FFT_T_N_i_Im;
                rx = FFT_T_i_Im   + FFT_T_N_i_Im;
                ry = FFT_T_N_i_Re - FFT_T_i_Re;

                lx_ = 2 * MagL[i] * Math.Cos(ArgL[i]);
                ly_ = 2 * MagL[i] * Math.Sin(ArgL[i]);
                rx_ = 2 * MagR[i] * Math.Cos(ArgR[i]);
                ry_ = 2 * MagR[i] * Math.Sin(ArgR[i]);

                currentDiff = Math.Abs(lx - lx_);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

                currentDiff = Math.Abs(ly - ly_);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

                currentDiff = Math.Abs(rx - rx_);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

                currentDiff = Math.Abs(ry - ry_);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
            }

            // Сохраняем максимальную невязку перехода из алгебраической формы
            // в показательную и обратно
            selfTestResult.MaxDiff_ALG_to_EXP_to_ALG = maxDiff;

            // Не используем взвешивающее окно, но работаем
            // с нормализацией - направление обратное
            useTaperWindow = false;
            FFT_S_Offset   = 0;
            recoverAfterTaperWindow = false;
            useNorm      = true;
            direction    = false;
            usePolyphase = false;
            CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
                         recoverAfterTaperWindow, useNorm, direction,
                         usePolyphase, fftObj);

            maxDiff = 0;
            for(i = 0; i < fftObj.N; ++i)
            {
                currentDiff = Math.Abs(FFT_S_backward[i] - FFT_S[i]);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
            }

            // Сохраняем максимальную невязку после прямого-обратного преобразования
            // (без взвешивающего окна)
            selfTestResult.MaxDiff_FORWARD_BACKWARD = maxDiff;

            // Используем взвешивающее окно, работаем
            // с нормализацией - направление прямое
            useTaperWindow = true;
            FFT_S_Offset   = 0;
            recoverAfterTaperWindow = false;
            useNorm      = true;
            direction    = true;
            usePolyphase = false;
            CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                         recoverAfterTaperWindow, useNorm, direction,
                         usePolyphase, fftObj);

            // Используем аннигиляцию взвешивающего окна, работаем
            // с нормализацией - направление обратное
            useTaperWindow = true;
            FFT_S_Offset   = 0;
            recoverAfterTaperWindow = true;
            useNorm      = true;
            direction    = false;
            usePolyphase = false;
            CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
                         recoverAfterTaperWindow, useNorm, direction,
                         usePolyphase, fftObj);

            maxDiff = 0;
            for(i = (fftObj.NN / 2); i <= ((fftObj.NN * 3) / 4); ++i)
            {
                currentDiff = Math.Abs(FFT_S_backward[i] - FFT_S[i]);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
            }

            // Сохраняем максимальную невязку после прямого-обратного
            // преобразования (c аннигиляцией взвешивающего окна)
            selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW = maxDiff;

            maxDiff = 0;
            for(i = 0; i < N2; ++i)
            {
                currentDiff = Math.Abs(ACH[i] - ACH_Difference);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
            }

            // Сохраняем максимальную невязку по расчету заданной АЧХ
            selfTestResult.MaxDiff_ACH = maxDiff;

            maxDiff = 0;
            for(i = 0; i < N2; ++i)
            {
                currentDiff = Math.Abs(PhaseLR[i]);
                maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
            }

            // Сохраняем максимальную невязку по расчету разности хода фаз каналов
            selfTestResult.MaxDiff_PhaseLR = maxDiff;

            // Performance Test
            timerFrequency = 10000000;
            CFFT_Process_counter = CFFT_Explore_counter = 0;

            // Не используем взвешивающее окно, но работаем
            // с нормализацией - направление прямое
            useTaperWindow = false;
            FFT_S_Offset   = 0;
            recoverAfterTaperWindow = false;
            useNorm      = false;
            direction    = true;
            usePolyphase = false;

            // CFFT_Process_time
            startCounter = 0;
            startCounter = DateTime.Now.Ticks;
            for(i = 0; i < N_iters; ++i)
            {
                CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                             recoverAfterTaperWindow, useNorm, direction,
                             usePolyphase, fftObj);
            }
            CFFT_Process_counter  = DateTime.Now.Ticks;
            CFFT_Process_counter -= startCounter;
            selfTestResult.CFFT_Process_time  = (double)CFFT_Process_counter / (double)timerFrequency;
            selfTestResult.CFFT_Process_time /= (double)N_iters;

            // CFFT_Explore_time
            startCounter = 0;
            startCounter = DateTime.Now.Ticks;
            for(i = 0; i < N_iters; ++i)
            {
                CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                             usePolyphase, fftObj);
            }
            CFFT_Explore_counter  = DateTime.Now.Ticks;
            CFFT_Explore_counter -= startCounter;
            selfTestResult.CFFT_Explore_time  = (double)CFFT_Explore_counter / (double)timerFrequency;
            selfTestResult.CFFT_Explore_time /= (double)N_iters;

            // Высвобождаем ресурсы динамической памяти
            FFT_S = null;
            FFT_S_backward = null;
            FFT_T   = null;
            MagL    = null;
            MagR    = null;
            ACH     = null;
            ArgL    = null;
            ArgR    = null;
            PhaseLR = null;

            // Проверка на допустимость полученных погрешностей
            if(selfTestResult.MaxDiff_ACH                     <= MAX_FFT_DIFF &&
               selfTestResult.MaxDiff_ALG_to_EXP_to_ALG       <= MAX_FFT_DIFF &&
               selfTestResult.MaxDiff_FORWARD_BACKWARD        <= MAX_FFT_DIFF &&
               selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW <= MAX_FFT_DIFF &&
               selfTestResult.MaxDiff_PhaseLR                 <= MAX_FFT_DIFF)
            {
                selfTestResult.AllOK = 1;
            }
            else
            {
                selfTestResult.AllOK = 0;
            }

            // Если активирован режим сброса дампа...
            if(IsDumpMode)
            {
                // Вердикт по точности выполнения самодиагностики
                DebugHelper.WriteInt(DumpName, "AllOK.int32", selfTestResult.AllOK);
                
                // Максимальная невязка по расчету заданной АЧХ
                DebugHelper.WriteDouble(DumpName, "MaxDiff_ACH.double", selfTestResult.MaxDiff_ACH);

                // Max. невязка ALG . EXP и обратно
                DebugHelper.WriteDouble(DumpName, "MaxDiff_ALG_to_EXP_to_ALG.double", selfTestResult.MaxDiff_ALG_to_EXP_to_ALG);

                // Max. невязка FORVARD + BACKWARD
                DebugHelper.WriteDouble(DumpName, "MaxDiff_FORWARD_BACKWARD.double", selfTestResult.MaxDiff_FORWARD_BACKWARD);

                //...то же + восст. после TW
                DebugHelper.WriteDouble(DumpName, "MaxDiff_FORWARD_BACKWARD_AntiTW.double", selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW);

                // Макс. невязка по расчету разности хода фаз
                DebugHelper.WriteDouble(DumpName, "MaxDiff_PhaseLR.double", selfTestResult.MaxDiff_PhaseLR);

                // Время работы CFFT_Process()
                DebugHelper.WriteDouble(DumpName, "CFFT_Process_time.double", selfTestResult.CFFT_Process_time);

                // Время работы CFFT_Explore()
                DebugHelper.WriteDouble(DumpName, "CFFT_Explore_time.double", selfTestResult.CFFT_Explore_time);
            }

            // Возвращаем результаты тестирования
            return selfTestResult;
        }

        /// <summary>
        /// Самотестирование внутренней точности в цикле прямого-обратного
        /// преобразования на случайных данных ("белый шум")
        /// </summary>
        /// <param name="ACH_Difference"> Коэффициент превосходства правого канала
        /// по уровню в ходе проводимого теста. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Структура "Результат тестирования внутренней точности
        /// прямого-обратного комплексного преобразования Фурье". </returns>
        public static CFFT_SelfTestResult SelfTest_RND(double ACH_Difference,
                                                       CFFT_Object fftObj)
        {
            int i;
            double randMult, randomValue;
            double[] FFT_S;
            Random rnd;

            // Массив исходных данных - стартовый
            FFT_S = new double[fftObj.NN];

            // Инициализируем генератор случайных чисел
            rnd = new Random();

            // Множитель случайной выборки
            randMult = 1E07;

            // Заполняем исходный массив случайными данными...
            for(i = 0; i < fftObj.N; ++i)
            {
                // Получаем значение с датчика случайных чисел
                //(эмулируем ввод белого шума)...
                randomValue = (rnd.NextDouble() * randMult) -
                              (rnd.NextDouble() * randMult);

                // "левый" канал и "правый" канал будут отличаться
                // в "ACH_Difference" раз
                FFT_S[(i << 1) + 0] = randomValue / ACH_Difference;
                FFT_S[(i << 1) + 1] = randomValue;
            }

            // Возвращаем результаты тестирования...
            return SelfTest_S(FFT_S, ACH_Difference, fftObj);
        }

        /// <summary>
        /// Поиск индекса максимального значения в массиве
        /// </summary>
        /// <param name="data"> Исходный массив для поиска максимума. </param>
        /// <param name="startIdx"> Частота семплирования. </param>
        /// <param name="finishIdx"> Глубина поиска. </param>
        public static int GetMaxIdx(double[] data, int startIdx, int finishIdx)
        {
            int i;
            int maxValIdxL, maxValIdxR;
            double currVal, maxVal;

            // Стартуем с первоначальным предположением...
            maxValIdxL = maxValIdxR = startIdx;
            maxVal     = data[maxValIdxL];

            for(i = startIdx + 1; i <= finishIdx; ++i)
            {
                currVal = data[i];

                // Если текущее больше максимума -
                // сдвигаем оба индекса вправо...
                if(currVal > maxVal)
                {
                    maxValIdxL = maxValIdxR = i;
                    maxVal     = currVal;
                }
                else
                {
                    //...а при равенстве максимуму
                    // сдвигается только правый индекс
                    if(currVal == maxVal)
                    {
                        maxValIdxR = i;
                    }
                }
            }

            return (maxValIdxL + maxValIdxR) >> 1;
        }

        /// <summary>
        /// Метод точного вычисления частоты по множеству гармоник
        /// </summary>
        /// <param name="Mag"> Магнитуды. </param>
        /// <param name="L"> Левая включенная граница для анализа. </param>
        /// <param name="R"> Правая включенная граница для анализа. </param>
        /// <param name="sampFreq"> Частота семплирования. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
        public static double CalcExactFreq(double[] Mag, int L, int R,
                                           double sampFreq,
                                           CFFT_Object fftObj)
        {
            int i;
            double harmSum, exactFreqIdx;

            // Узнаем сумму гармоник
            harmSum = 0;
            for(i = L; i <= R; ++i)
            {
                harmSum += Mag[i];
            }

            // Вычисляем индекс точной частоты
            exactFreqIdx = 0;
            for(i = L; i <= R; ++i)
            {
                exactFreqIdx += (Mag[i] / harmSum) * (double)i;
            }

            // Возвращаем точную частоту
            return FreqNode(exactFreqIdx, sampFreq, fftObj.N, fftObj.IsComplex);
        }

        /// <summary>
        /// Метод точного вычисления частоты по множеству гармоник
        /// </summary>
        /// <param name="Mag"> Магнитуды. </param>
        /// <param name="L"> Левая включенная граница для анализа. </param>
        /// <param name="R"> Правая включенная граница для анализа. </param>
        /// <param name="depth"> Глубина поиска. </param>
        /// <param name="sampFreq"> Частота семплирования. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
        public static double ExactFreq(double[] Mag, int L, int R, int depth,
                                       double sampFreq,
                                       CFFT_Object fftObj)
        {
            int fft_N, FFT_NodeMax, startIdx, finishIdx, startIdx2, finishIdx2,
                deltaStart, deltaFinish, deltaCorr;

            // Ищем максимальный индекс в установленных границах
            FFT_NodeMax = GetMaxIdx(Mag, L, R);

            // Вычисляем первоначальные индексы
            startIdx  = L;
            finishIdx = R;

            // Количество гармоник FFT
            fft_N = fftObj.N * (fftObj.IsComplex ? 2 : 1);

            // Корректируем индексы
            startIdx2  = startIdx  >= 0 ? startIdx : 0;
            finishIdx2 = finishIdx <= (fft_N - 1) ? finishIdx : (fft_N - 1);

            // Глубину обработки корректируем на максимальную дельту индексов
            deltaStart  = Math.Abs(startIdx2  - startIdx);
            deltaFinish = Math.Abs(finishIdx2 - finishIdx);
            deltaCorr   = Math.Max(deltaStart, deltaFinish);
            depth -= deltaCorr;

            // Вычисляем окончательные индексы
            startIdx  = FFT_NodeMax - depth;
            finishIdx = FFT_NodeMax + depth;

            // Возвращаем точную частоту
            return CalcExactFreq(Mag, startIdx, finishIdx, sampFreq, fftObj);
        }

        /// <summary>
        /// Метод точного вычисления частоты по множеству гармоник
        /// </summary>
        /// <param name="Mag">  Магнитуды. </param>
        /// <param name="depth"> Глубина поиска. </param>
        /// <param name="sampFreq"> Частота семплирования. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
        public static double ExactFreqAuto(double[] Mag, int depth,
                                           double sampFreq,
                                           CFFT_Object fftObj)
        {
            return ExactFreq(Mag, 1, ((fftObj.N * (fftObj.IsComplex ? 2 : 1)) >> 1) - 1, depth, sampFreq, fftObj);
        }

        /// <summary>
        /// Преобразование вещественного входа в комплексный
        /// </summary>
        /// <param name="Re"> Вещественная компонента. </param>
        /// <param name="Im"> Мнимая компонента. </param>
        /// <returns> Комплексный выход. </returns>
        public static double[] RealToComplex(double[] Re, double[] Im)
        {
            int i, j;
            double[] Complex;

            if(Re.Length != Im.Length)
            {
                throw new Exception("ExactFFT::RealToComplex(): (Re.Length != Im.Length)");
            }

            Complex = new double[Re.Length << 1];
            for(i = 0, j = 0; i < Re.Length; i++, j += 2)
            {
                Complex[j + 0] = Re[i];
                Complex[j + 1] = Im[i];
            }

            return Complex;
        }

        /// <summary>
        /// Преобразование вещественного входа в комплексный
        /// </summary>
        /// <param name="real"> Вещественный вход. </param>
        /// <param name="complexPart"> Выбор целевой части преобразования:
        ///                            real <=> 0, img <=> 1. </param>
        /// <returns> Комплексный выход. </returns>
        public static double[] RealToComplex(this double[] real, int complexPart)
        {
            int i, j;
            var complex = new double[real.Length << 1];
            for(i = 0, j = 0; i < real.Length; i++, j += 2)
            {
                complex[j + ((complexPart + 0) % 2)] = real[i];
                complex[j + ((complexPart + 1) % 2)] = 0;
            }

            return complex;
        }

        /// <summary>
        /// Преобразование комплексного входа в вещественный
        /// </summary>
        /// <param name="complex"> Комплексный вход. </param>
        /// <param name="complexPart"> Выбор целевой части преобразования:
        ///                            real <=> 0, img <=> 1. </param>
        /// <returns> Вещественный выход. </returns>
        public static double[] ComplexToReal(this double[] complex, int complexPart)
        {
            int i, j;
            var real = new double[complex.Length >> 1];
            for(i = 0, j = 0; i < real.Length; i++, j += 2)
            {                
                real[i] = complex[j + (complexPart % 2)];
            }

            return real;
        }                
    }
}