/*----------------------------------------------------------------------+
 |  filename:   ExactPlotter.cs                                         |
 |----------------------------------------------------------------------|
 |  version:    8.30                                                    |
 |  revision:   13/01/2016  21:27                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Многопоточный "плоттер-водопад" для комплексного FFT    |
 |----------------------------------------------------------------------*/

using System;
using System.Threading.Tasks;

namespace DrAF.DSP
{
    /// <summary>
    /// Многопоточный "плоттер-водопад" для комплексного FFT
    /// </summary>
    public class ExactPlotter
    {
        /// <summary>
        /// Результат разбора выходных данных CFFT
        /// </summary>
        public class CFFT_ExploreResult
        {
            //--- Explore() -----------------------------------------------------------
            public double[][] MagL;    // Магнитуды "левого" канала.
            public double[][] MagR;    // Магнитуды "правого" канала.
            public double[][] ACH;     // АЧХ (отношение магнитуды "правого" канала к магнитуде "левого" - как "выход" / "вход").
            public double[][] ArgL;    // Аргумент "левого" канала.
            public double[][] ArgR;    // Аргумент "правого" канала.
            public double[][] PhaseLR; // Разность хода фаз каналов ("правый" минус "левый").
            //--- ComplexExplore() ----------------------------------------------------
            public double[][] Mag;     // Магнитуды.
            public double[][] Arg;     // Аргументы.
        }

        /// <summary>
        /// Вычисление количества строк в "водопаде"
        /// </summary>
        /// <param name="FFT_S"> Вектор входных данных
        /// ("левый" и "правый" каналы - чет./нечет.). </param>
        /// <param name="FFT_S_Offset"> Смещение данных для анализа во
        /// входном векторе FFT_S. </param> 
        /// <param name="remainArrayItemsLRCount"> Остаток необработанных данных в исходном
        /// массиве (количество элементов, включая Re/Im). </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static int GetPlotRowsCount(int FFT_S_Length, int FFT_S_Offset,
                                           out int remainArrayItemsLRCount,
                                           ExactFFT.CFFT_Object fftObj)
        {
            int stepAreaLength_Real, nSteps, areaWithoutFFT_Real, intersectFFTArea_Real;

            // Вычисляем размер поля данных для "поглощения" шагами FFT.
            // Деление на 2 требуется для перехода от комплексной "чет-нечет" формы к "вещественной"
            stepAreaLength_Real = ((FFT_S_Length - FFT_S_Offset) >> 1) - fftObj.N;

            // Невозможно обработать столь малый блок (некуда шагать)!
            if(stepAreaLength_Real < 0)
            {
                throw new Exception("ExactPlotter::GetPlotRowsCount: (stepAreaLength_Real < 0)");
            }

            // Количество шагов взвешивающим окном равно размеру области, первоначально не лежащей
            // под окном FFT, деленному на шаг окна FFT. Так-как шаги делаются целиком, нельзя
            // осуществить полшага, поэтому берется только целая часть числа.
            // Оставшееся количество попадает в необработанный остаток!
            // +1 учитывает "нулевой" шаг, при первоначальном положении на исходных данных.
            nSteps = ((stepAreaLength_Real / fftObj.WindowStep) + 1);

            // Когда вычисляется необработанный остаток, следует понимать, что следующий шаг FFT
            // покроет остаток, и выйдет даже за пределы поля исходных данных, поэтому процесс
            // и останавливается. Но если сделать шаг FFT, то оно захватит не только всё поле
            // необработанных данных, если его шаг меньше его самого, произойдет частичный
            // захват области под ним самим (в состоянии итерации при останове).
            // Размер этой области равен размеру окна FFT минус его шаг. И эта величина также
            // требует, чтобы её учитывали!
            // Если шаг FFT больше его самого, получится отрицательная величина.
            areaWithoutFFT_Real = stepAreaLength_Real - ((nSteps - 1) * fftObj.WindowStep);
            intersectFFTArea_Real = fftObj.N - fftObj.WindowStep;
            remainArrayItemsLRCount = (intersectFFTArea_Real + areaWithoutFFT_Real) << 1;

            // Отрицательный остаток обработки дает ошибку при построении сонограммы!
            if(remainArrayItemsLRCount < 0)
            {
                throw new Exception("ExactPlotter::GetPlotRowsCount: (remainArrayItemsLRCount < 0)");
            }

            return nSteps;
        }


        /// <summary>
        /// Вычисление количества строк в "водопаде"
        /// </summary>
        /// <param name="FFT_S"> Вектор входных данных
        /// ("левый" и "правый" каналы - чет./нечет.). </param>
        /// <param name="FFT_S_Offset"> Смещение данных для анализа во
        /// входном векторе FFT_S. </param>
        /// <param name="remainArrayItemsLRCount"> Остаток необработанных данных в исходном
        /// массиве (количество элементов, включая Re/Im). </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static int GetPlotRowsCount(double[] FFT_S, int FFT_S_Offset,
                                           out int remainArrayItemsLRCount,
                                           ExactFFT.CFFT_Object fftObj)
        {
            if(FFT_S == null)
            {
                throw new Exception("ExactPlotter::GetPlotRowsCount(): (FFT_S == null)");
            }

            return GetPlotRowsCount(FFT_S.Length, FFT_S_Offset, out remainArrayItemsLRCount, fftObj);
        }
        
        /// <summary>
        /// Построение комплексной сонограммы
        /// </summary>
        /// <param name="FFT_S"> Вектор входных данных
        /// ("левый" и "правый" каналы - чет./нечет.). </param>
        /// <param name="FFT_S_Offset"> Смещение данных для анализа во
        /// входном векторе FFT_S. </param>
        /// <param name="useTaperWindow"> Использовать взвешивающее окно? </param>
        /// <param name="recoverAfterTaperWindow"> Аннигилировать действие
        /// взвешивающего окна на обратном проходе? </param>
        /// <param name="useNorm"> Использовать нормализацию 1/N? </param>
        /// <param name="direction"> Направление преобразования (true - прямое). </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="remainArrayItemsLRCount"> Остаток необработанных данных в исходном
        /// массиве (количество элементов, включая Re/Im). </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Сонограмма. </returns>
        public static double[][] Process(double[] FFT_S, int FFT_S_Offset,
                                         bool useTaperWindow, bool recoverAfterTaperWindow,
                                         bool useNorm, bool direction, bool usePolyphase,
                                         out int remainArrayItemsLRCount,
                                         ExactFFT.CFFT_Object fftObj)
        {
            int plotRowsCount, frameOffset;
            double[][] FFT_T;

            if(FFT_S == null)
            {
                throw new Exception("ExactPlotter::Process(): (FFT_S == null)");
            }

            // Вычисляем количество строк сонограммы...
            plotRowsCount = GetPlotRowsCount(FFT_S, FFT_S_Offset, out remainArrayItemsLRCount, fftObj);

            // Обрабатываем все фреймы...
            FFT_T = new double[plotRowsCount][];
            Parallel.For(0, plotRowsCount, frame =>
            {
                // Умножение на 2 треб. для real -> complex
                FFT_T[frame] = new double[fftObj.N << 1];
                frameOffset = FFT_S_Offset + frame * fftObj.WindowStep << 1;
                ExactFFT.CFFT_Process(FFT_S, frameOffset, FFT_T[frame], useTaperWindow,
                                      recoverAfterTaperWindow, useNorm, direction, usePolyphase,
                                      fftObj);
            });

            return FFT_T;
        }

        /// <summary>
        /// Извлечение данных из комплексной сонограммы (L+R)
        /// </summary>
        /// <param name="FFT_T"> Выходной набор векторов коэффициентов. </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Результат разбора выходных данных CFFT. </returns>
        public static CFFT_ExploreResult Explore(double[][] FFT_T, bool usePolyphase,
                                                 ExactFFT.CFFT_Object fftObj)
        {
            int plotRowsCount;
            CFFT_ExploreResult res;
            res = new CFFT_ExploreResult();

            if(FFT_T == null)
            {                
                throw new Exception("ExactPlotter::Explore(): (FFT_T == null)");
            }
            
            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = FFT_T.Length;

            // Подготавливаем выходные массивы...
            res.MagL    = new double[plotRowsCount][];
            res.MagR    = new double[plotRowsCount][];
            res.ACH     = new double[plotRowsCount][];
            res.ArgL    = new double[plotRowsCount][];
            res.ArgR    = new double[plotRowsCount][];
            res.PhaseLR = new double[plotRowsCount][];

            // Работаем по всем строкам сонограммы...
            Parallel.For(0, plotRowsCount, frame =>            
            {
                // Количество гармоник в два раза меньше размера кадра FFT
                res.MagL[frame]    = new double[fftObj.N >> 1];
                res.MagR[frame]    = new double[fftObj.N >> 1];
                res.ACH[frame]     = new double[fftObj.N >> 1];
                res.ArgL[frame]    = new double[fftObj.N >> 1];
                res.ArgR[frame]    = new double[fftObj.N >> 1];
                res.PhaseLR[frame] = new double[fftObj.N >> 1];
                
                // Извлечение данных FFT из комплексной сонограммы
                ExactFFT.CFFT_Explore(FFT_T[frame],
                                      res.MagL[frame], res.MagR[frame], res.ACH[frame],
                                      res.ArgL[frame], res.ArgR[frame], res.PhaseLR[frame],
                                      usePolyphase,
                                      fftObj);
            });

            return res;
        }

        /// <summary>
        /// Исследование результатов комплексного FFT (идентично CFFT из MathCAD)
        /// </summary>
        /// <param name="FFT_T"> Выходной набор векторов коэффициентов. </param>        
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="isMirror"> Зеркальное отображение спектра? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Результат разбора выходных данных CFFT. </returns>
        public static CFFT_ExploreResult ComplexExplore(double[][] FFT_T,
                                                        bool usePolyphase, bool isMirror,
                                                        ExactFFT.CFFT_Object fftObj)
        {
            int plotRowsCount;
            CFFT_ExploreResult res;
            res = new CFFT_ExploreResult();

            if(FFT_T == null)
            {
                throw new Exception("ExactPlotter::ComplexExplore(): (FFT_T == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = FFT_T.Length;

            // Подготавливаем выходные массивы...
            res.Mag = new double[plotRowsCount][];
            res.Arg = new double[plotRowsCount][];
            
            // Работаем по всем строкам сонограммы...
            Parallel.For(0, plotRowsCount, frame =>
            {
                // Количество гармоник равно размеру кадра FFT
                res.Mag[frame] = new double[fftObj.N];
                res.Arg[frame] = new double[fftObj.N];

                // Извлечение данных FFT из комплексной сонограммы (режим COMPLEX)
                ExactFFT.CFFT_ComplexExplore(FFT_T[frame],
                                             res.Mag[frame], res.Arg[frame],
                                             usePolyphase,
                                             isMirror,
                                             fftObj);
            });

            return res;
        }

        /// <summary>
        /// Исследование результатов комплексного FFT (идентично CFFT из MathCAD)
        /// </summary>
        /// <param name="FFT_T"> Выходной набор векторов коэффициентов. </param>
        /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
        /// <param name="isMirror"> Зеркальное отображение спектра? </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <returns> Результат разбора выходных данных CFFT. </returns>
        public static CFFT_ExploreResult ComplexExploreMag(double[][] FFT_T,
                                                           bool usePolyphase, bool isMirror,
                                                           ExactFFT.CFFT_Object fftObj)
        {
            int plotRowsCount;
            CFFT_ExploreResult res;
            res = new CFFT_ExploreResult();

            if(FFT_T == null)
            {
                throw new Exception("ExactPlotter::ComplexExploreMag(): (FFT_T == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = FFT_T.Length;

            // Подготавливаем выходные массивы...
            res.Mag = new double[plotRowsCount][];
            
            // Работаем по всем строкам сонограммы...
            Parallel.For(0, plotRowsCount, frame =>
            {
                // Количество гармоник равно размеру кадра FFT
                res.Mag[frame] = new double[fftObj.N];

                // Извлечение данных FFT из комплексной сонограммы (режим COMPLEX)
                ExactFFT.CFFT_ComplexExplore(FFT_T[frame],
                                             res.Mag[frame], null,
                                             usePolyphase,
                                             isMirror,
                                             fftObj);
            });

            return res;
        }

        /// <summary>
        /// Перевод значений массива double в форму dB
        /// </summary>
        /// <param name="sonogram"> Данные для обработки. </param>
        /// <param name="zero_db_level"> Значение "нулевого" уровня. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        public static void dB_Scale(double[][] sonogram, double zero_db_level,
                                    ExactFFT.CFFT_Object fftObj)
        {            
            int plotRowsCount;

            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::dB_Scale(): (sonogram == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = sonogram.Length;

            // Работаем по всем строкам сонограммы...
            Parallel.For(0, plotRowsCount, frame =>
            {
                ExactFFT.dB_Scale(sonogram[frame], zero_db_level, fftObj);
            });
        }

        /// <summary>
        /// Выделение одной гармоники из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="harmFreq"> Частота выделяемой гармоники. </param>
        /// <param name="harmIdx"> Индекс выделяемой гармоники. </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <param name="sampleRate"> Частота семплирования. </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[] HarmSlice(double[][] sonogram, double harmFreq, out int harmIdx,
                                         ExactFFT.CFFT_Object fftObj,
                                         int sampleRate)
        {
            double[] harmVector;
            double[][] harmColumn;

            // Выделяем ОДНУ гармонику (она размещается в одном столбце и множестве строк)...
            harmColumn = SubBand(sonogram, harmFreq, harmFreq, out harmIdx, out harmIdx, fftObj, sampleRate);

            // Данные вектора-столбца размещаем в векторе-строке...
            harmVector = new double[harmColumn.Length];

            // Работаем по всем строкам сонограммы...
            Parallel.For(0, harmColumn.Length, frame =>
            {
                harmVector[frame] = harmColumn[frame][0];
            });

            //...и возвращаем результат
            return harmVector;
        }

        /// <summary>
        /// Выделение одной гармоники из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="harmIdx"> Индекс выделяемой гармоники. </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[] HarmSlice(double[][] sonogram, int harmIdx)
        {
            double[] harmVector;

            // Выделяем ОДНУ гармонику (она размещается в одном столбце и множестве строк)...
            var harmColumn = SubBand(sonogram, harmIdx, harmIdx);

            // Данные вектора-столбца размещаем в векторе-строке...
            harmVector = new double[harmColumn.Length];

            // Работаем по всем строкам вектора-столбца...
            Parallel.For(0, harmColumn.Length, frame =>
            {
                harmVector[frame] = harmColumn[frame][0];
            });

            //...и возвращаем результат
            return harmVector;
        }

        /// <summary>
        /// Выделение поддиапазона гармоник из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="lowFreq"> Частота нижней гармоники выделяемого диапазона. </param>
        /// <param name="highFreq"> Частота верхней гармоники выделяемого диапазона. </param>
        /// <param name="lowHarmIdx"> Нижний индекс гармоники. </param>
        /// <param name="highHarmIdx"> Верхний индекс гармоники. </param>        
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <param name="sampleRate"> Частота семплирования. </param>
        /// <param name="harmReverse"> Требуется реверс гармоник? </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[][] SubBand(double[][] sonogram, double lowFreq, double highFreq,
                                         out int lowHarmIdx, out int highHarmIdx,
                                         ExactFFT.CFFT_Object fftObj,
                                         int sampleRate,
                                         bool harmReverse = false)
        {
            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::SubBand(): (sonogram == null)");
            }

            lowHarmIdx  = (int)ExactFFT.FFT_Node(lowFreq,  sampleRate, fftObj.N, fftObj.IsComplex);
            highHarmIdx = (int)ExactFFT.FFT_Node(highFreq, sampleRate, fftObj.N, fftObj.IsComplex);
            
            return SubBand(sonogram, lowHarmIdx, highHarmIdx, harmReverse);
        }

        /// <summary>
        /// Выделение поддиапазона гармоник из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="lowHarmIdx"> Нижний индекс гармоники. </param>
        /// <param name="highHarmIdx"> Верхний индекс гармоники. </param>
        /// <param name="harmReverse"> Требуется реверс гармоник? </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[][] SubBand(double[][] sonogram, int lowHarmIdx, int highHarmIdx,
                                         bool harmReverse = false)
        {
            int plotRowsCount, plotColsCount;
            double[][] target;

            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::SubBand(): (sonogram == null)");
            }

            // Считываем количество столбцов, которые имеет сонограмма...
            plotColsCount = sonogram[0].Length;

            // Поддиапазон должен состоять как минимум из одного элемента
            if((highHarmIdx - lowHarmIdx) < 0)
            {
                throw new Exception(("ExactPlotter::SubBand(): ((highHarmIdx - lowHarmIdx) < 0)"));
            }

            // Проверка на допустимость индекса "startFrameIdx"
            if((lowHarmIdx < 0) || (lowHarmIdx >= plotColsCount))
            {
                throw new Exception(("ExactPlotter::SubBand(): ((lowHarmIdx < 0) || (lowHarmIdx >= plotColsCount))"));
            }

            // Проверка на допустимость индекса "finishFrameIdx"
            if((highHarmIdx < 0) || (highHarmIdx >= plotColsCount))
            {
                throw new Exception(("ExactPlotter::SubBand(): ((highHarmIdx < 0) || (highHarmIdx >= plotColsCount))"));
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = sonogram.Length;

            // Целевое хранилище сонограммы
            target = new double[plotRowsCount][];
            
            // Работаем по всем строкам сонограммы...
            Parallel.For(0, plotRowsCount, frame =>
            {
                double[] sourceRow = sonogram[frame];
                double[] targetRow = new double[(highHarmIdx - lowHarmIdx) + 1];

                Array.Copy(sourceRow, lowHarmIdx, targetRow, 0, (highHarmIdx - lowHarmIdx) + 1);
                if(harmReverse)
                {
                    Array.Reverse(targetRow);
                }

                target[frame] = targetRow;
            });

            return target;
        }

        /// <summary>
        /// Выделение одного момента времени из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="time"> Время среза (кадра FFT). </param>
        /// <param name="timeIdx"> Индекс временного среза (кадра FFT). </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <param name="sampleRate"> Частота семплирования. </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[] TimeSlice(double[][] sonogram, double time, out int timeIdx,
                                         ExactFFT.CFFT_Object fftObj,
                                         int sampleRate)
        {
            // Выделяем один временной срез (выделяем вектор, привязанный к ОДНОЙ строке)...
            return SubTime(sonogram, time, time, out timeIdx, out timeIdx, fftObj, sampleRate)[0];
        }

        /// <summary>
        /// Выделение одного момента времени из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="timeIdx"> Индекс временного среза (кадра FFT). </param>        
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[] TimeSlice(double[][] sonogram, int timeIdx)
        {
            // Выделяем один временной срез (выделяем вектор, привязанный к ОДНОЙ строке)...
            return SubTime(sonogram, timeIdx, timeIdx)[0];
        }
        
        /// <summary>
        /// Выделение поддиапазона времени из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="startTime"> Стартовый кадр FFT (время строки сонограммы). </param>
        /// <param name="finishTime"> Конечный кадр FFT (время строки сонограммы). </param>
        /// <param name="startFrameIdx"> Стартовый кадр FFT (строка сонограммы). </param>
        /// <param name="finishFrameIdx"> Конечный кадр FFT (строка сонограммы). </param>
        /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
        /// <param name="sampleRate"> Частота семплирования. </param>
        /// <param name="harmReverse"> Требуется реверс гармоник? </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[][] SubTime(double[][] sonogram, double startTime, double finishTime,
                                         out int startFrameIdx, out int finishFrameIdx,
                                         ExactFFT.CFFT_Object fftObj,
                                         int sampleRate,
                                         bool harmReverse = false)
        {
            double frameDuration, sleepCoeff, timeSliceDuration;
            GetFrameParameters(fftObj.N, fftObj.WindowStep, sampleRate,
                               out frameDuration, out sleepCoeff, out timeSliceDuration);

            startFrameIdx  = (int)(startTime  / timeSliceDuration);
            finishFrameIdx = (int)(finishTime / timeSliceDuration);

            return SubTime(sonogram, startFrameIdx, finishFrameIdx, harmReverse);
        }

        /// <summary>
        /// Выделение поддиапазона времени из сонограммы
        /// </summary>
        /// <param name="sonogram"> Исходная сонограмма. </param>
        /// <param name="startFrameIdx"> Стартовый кадр FFT (строка сонограммы). </param>
        /// <param name="finishFrameIdx"> Конечный кадр FFT (строка сонограммы). </param>
        /// <param name="harmReverse"> Требуется реверс гармоник? </param>
        /// <returns> Выделенный поддиапазон сонограммы. </returns>
        public static double[][] SubTime(double[][] sonogram, int startFrameIdx, int finishFrameIdx,
                                         bool harmReverse = false)
        {            
            int plotRowsCount;
            double[][] target;

            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::SubTime(): (sonogram == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = sonogram.Length;

            // В выделяемом фрагменте времени должен быть по крайней мере один фрейм
            if((finishFrameIdx - startFrameIdx) < 0)
            {
                throw new Exception(("ExactPlotter::SubTime(): ((finishFrameIdx - startFrameIdx) < 0)"));
            }

            // Проверка на допустимость индекса "startFrameIdx"
            if((startFrameIdx < 0) || (startFrameIdx >= plotRowsCount))
            {
                throw new Exception(("ExactPlotter::SubTime(): ((startFrameIdx < 0) || (startFrameIdx >= plotRowsCount))"));
            }

            // Проверка на допустимость индекса "finishFrameIdx"
            if((finishFrameIdx < 0) || (finishFrameIdx >= plotRowsCount))
            {
                throw new Exception(("ExactPlotter::SubTime(): ((finishFrameIdx < 0) || (finishFrameIdx >= plotRowsCount))"));
            }
            
            // Целевое хранилище сонограммы
            target = new double[(finishFrameIdx - startFrameIdx) + 1][];

            // Работаем по всем строкам сонограммы...            
            Parallel.For(startFrameIdx, (finishFrameIdx + 1), frame =>
            {
                double[] sourceRow = sonogram[frame];
                double[] targetRow = new double[sourceRow.Length];
                Array.Copy(sourceRow, targetRow, sourceRow.Length);

                if(harmReverse)
                {
                    Array.Reverse(targetRow);
                }

                target[frame - startFrameIdx] = targetRow;
            });

            return target;
        }

        /// <summary>
        /// Получение временнЫх параметров FFT
        /// </summary>
        /// <param name="windowSize"> Размер окна FFT. </param>
        /// <param name="windowStep"> Шаг окна FFT. </param>
        /// <param name="sampleRate"> Частота семплирования. </param>
        /// <param name="frameDuration"> Длительность кадра FFT. </param>
        /// <param name="sleepCoeff"> Коэффициент "скольжения". </param>
        /// <param name="timeSliceDuration"> Размер "временнОго среза". </param>
        public static void GetFrameParameters(int windowSize, int windowStep, int sampleRate,
                                              out double frameDuration, out double sleepCoeff,
                                              out double timeSliceDuration)
        {
            frameDuration     = (double)windowSize / sampleRate; // Длительность кадра FFT
            sleepCoeff        = (double)windowSize / windowStep; // Коэффициент "скольжения"
            timeSliceDuration = frameDuration      / sleepCoeff; // Длительность "временного среза"
        }

        /// <summary>
        /// Почленная сумма двух векторов, размещается в первом векторе
        /// </summary>
        /// <param name="summa"> Вектор-аккумулятор. </param>
        /// <param name="vector"> Вектор-член, входящий в сумму. </param>
        public static void VectorSum(double[] summa, double[] vector)
        {
            int i;

            if(summa == null)
            {
                throw new Exception("ExactPlotter::VectorSum(): (summa == null)");
            }

            if(vector == null)
            {
                throw new Exception("ExactPlotter::VectorSum(): (vector == null)");
            }

            if(summa.Length != vector.Length)
            {
                throw new Exception("ExactPlotter::VectorSum(): (summa.Length != vector.Length)");
            }

            for(i = 0; i < vector.Length; i++)
            {
                summa[i] += vector[i];
            }
        }
        
        /// <summary>
        /// Вычисление средневзвешенной проекции сонограммы на ось частот
        /// </summary>
        /// <param name="sonogram"> Сонограмма или её выделенный фрагмент. </param>
        /// <returns> Средневзвешенная проекция сонограммы на ось частот. </returns>
        public static double[] TimeSlicesSum(double[][] sonogram)
        {
            int plotRowsCount, plotColsCount;
            double div;
            double[] timeSlice, summa;

            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::TimeSlicesSum(): (sonogram == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = sonogram.Length;

            // Считываем количество столбцов, которые имеет сонограмма...
            plotColsCount = sonogram[0].Length;

            // В векторе суммы столько элементов, сколько гармоник в сонограмме...
            summa = new double[plotColsCount];
            for(int frame = 0; frame < plotRowsCount; frame++)
            {
                timeSlice = TimeSlice(sonogram, frame);
                VectorSum(summa, timeSlice);
            }

            // Нормализуем результат (по количеству слагаемых в векторе)...
            div = (double)plotRowsCount;
            Parallel.For(0, summa.Length, i =>
            {
                summa[i] /= div;
            });
           
            return summa;
        }

        /// <summary>
        /// Вычисление средневзвешенной проекции сонограммы на ось времени
        /// </summary>
        /// <param name="sonogram"> Сонограмма или её выделенный фрагмент. </param>
        /// <returns> Средневзвешенная проекция сонограммы на ось времени. </returns>
        public static double[] HarmSlicesSum(double[][] sonogram)
        {
            int plotRowsCount, plotColsCount;
            double div;
            double[] harmSlice, summa;

            if(sonogram == null)
            {
                throw new Exception("ExactPlotter::HarmSlicesSum(): (sonogram == null)");
            }

            // Считываем количество строк, которые имеет сонограмма...
            plotRowsCount = sonogram.Length;

            // Считываем количество столбцов, которые имеет сонограмма...
            plotColsCount = sonogram[0].Length;
            
            // В векторе суммы столько элементов, сколько фреймов в сонограмме...
            summa = new double[plotRowsCount];
            for(int harm = 0; harm < plotColsCount; harm++)
            {
                harmSlice = HarmSlice(sonogram, harm);
                VectorSum(summa, harmSlice);
            }

            // Нормализуем результат (по количеству слагаемых в векторе)...
            div = (double)plotColsCount;
            Parallel.For(0, summa.Length, i =>
            {
                summa[i] /= div;
            });

            return summa;
        }        
    }
}
