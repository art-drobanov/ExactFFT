/*----------------------------------------------------------------------+
 |  filename:   ExactFFT.c                                              |
 |----------------------------------------------------------------------|
 |  version:    7.10                                                    |
 |  revision:   07/09/2014  11:41                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Комплексное FFT                                         |
 |----------------------------------------------------------------------*/

 #ifndef _exactfft_c
 #define _exactfft_c

 #include "ExactFFT.h"
 #include "Windows.h"

 extern enum CosTW;

 /// <summary>
 /// Получение указателя на файловый поток
 /// <summary>
 /// <param name="dirName"> Имя директории. </param>
 /// <param name="fileName"> Имя файла. </param>
 /// <param name="toWrite"> Открыть файл для записи? </param>
 FILE * GetStreamPointer(char *dirName, char *fileName, bool toWrite)
 {
     // Полный путь к файлу дампа
     char path[MAX_PATH];

     // Подготавливаем полный путь
     if(dirName != NULL)
     {
        strcpy(path, dirName);
        strcat(path, "\\");
        strcat(path, fileName);

     } else
     {
        strcpy(path, fileName);
     }

     // Работа с файлом
     return fopen(path, toWrite ? "wb" : "rb");
 }

 /// <summary>
 /// Сброс массива int-ов в файл в двоичной форме
 /// <summary>
 /// <param name="arr"> Исходный массив. </param>
 /// <param name="N"> Количество элементов для сброса. </param>
 /// <param name="dirName"> Имя директории дампа. </param>
 /// <param name="fileName"> Имя файла дампа. </param>
 void DumpInt(int *arr, int N, char *dirName, char *fileName)
 {
     int i;
     int data;
     FILE *f = GetStreamPointer(dirName, fileName, TRUE);

     for(i = 0; i < N; ++i)
     {
         data = arr[i];
         fwrite(&data, sizeof(int), 1, f);
     }

     fflush(f); fclose(f);
 }

 /// <summary>
 /// Сброс массива double-ов в файл в двоичной форме
 /// <summary>
 /// <param name="arr"> Исходный массив. </param>
 /// <param name="N"> Количество элементов для сброса. </param>
 /// <param name="dirName"> Имя директории дампа. </param>
 /// <param name="fileName"> Имя файла дампа. </param>
 void DumpDouble(double *arr, int N, char *dirName, char *fileName)
 {
     int i;
     double data;
     FILE *f = GetStreamPointer(dirName, fileName, TRUE);

     for(i = 0; i < N; ++i)
     {
         data = arr[i];
         fwrite(&data, sizeof(double), 1, f);
     }

     fflush(f); fclose(f);
 }

 /// <summary>
 /// Логарифм по произвольному основанию
 /// </summary>
 /// <param name="arg"> Аргумент логарифма. </param>
 /// <param name="logBase"> Основание логарифма. </param>
 double LogX(double arg, double logBase)
 {
     return log(arg) / log(logBase);
 }

 /// <summary>
 /// Приведение значения к ближайшей снизу степени двойки
 /// </summary>
 /// <param name="arg"> Входной аргумент. </param>
 int ToLowerPowerOf2(int arg)
 {
     return (int)pow(2, (int)(LogX(arg, 2)));
 }

 /// <summary>
 /// Получение частоты заданного узла FFT
 /// </summary>
 /// <param name="FFT_Node"> Номер гармоники. </param>
 /// <param name="sampFreq"> Частота семплирования. </param>
 /// <param name="isComplex"> Комплексный режим? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 double FreqNode(double FFT_Node, double sampFreq, bool isComplex,
                 CFFT_Object *fftObj)
 {
     return (FFT_Node * sampFreq) / ((double)fftObj->N * (isComplex ? 2.0 : 1.0));
 }

 /// <summary>
 /// Получение узла FFT по заданной частоте
 /// </summary>
 /// <param name="freqNode"> Заданная частота. </param>
 /// <param name="sampFreq"> Частота семплирования. </param>
 /// <param name="isComplex"> Комплексный режим? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 double FFT_Node(double freqNode, double sampFreq, bool isComplex,
                 CFFT_Object *fftObj)
 {
     return (freqNode * ((double)fftObj->N * (isComplex ? 2.0 : 1.0))) / sampFreq;
 }

 /// <summary>
 /// Нормирование фазы
 /// </summary>
 /// <param name="phase"> Значение фазы для нормирования. </param>
 double PhaseNorm(double phase)
 {
     if(phase > 0) while (phase >= M_PI)  phase -= M_2PI;
              else while (phase <= -M_PI) phase += M_2PI;
 
     return phase;
 }

 /// <summary>
 /// Безопасное в смысле значений аргументов вычисление арктангенса
 /// </summary>
 /// <param name="im"> Мнимая часть комплексного числа. </param>
 /// <param name="re"> Действительная часть комплексного числа. </param>
 double Safe_atan2(double im, double re)
 {
     return (((re < 0) ? -re : re) < FLOAT_MIN) ? 0 : atan2(im, re);
 }

 /// <summary>
 /// Заполнение вектора изменения порядка следования данных перед FFT
 /// </summary>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void fill_FFT_P(CFFT_Object *fftObj)
 {
     int i, j, shift;

     // Выделяем память под вектор перестановок FFT...
     SAFE_DELETE(fftObj->FFT_P);
     fftObj->FFT_P = (int *)calloc(fftObj->NN, sizeof(int));

     // Инициализация массива перестановок...
     memset(fftObj->FFT_P, 0x00, (fftObj->NN * sizeof(int)));

     // Заполняем вектор изменения порядка следования данных...
     for(j = 0; j < LogX(fftObj->N, 2); ++j)
     {
        for(i = 0; i < fftObj->N; ++i)
        {
            fftObj->FFT_P[i << 1] = ((fftObj->FFT_P[i << 1] << 1) +
                                    ((i >> j) & 1));
        }
     }

     shift = (fftObj->FFT_P[2] == (fftObj->N >> 1)) ? 1 : 0;
     for(i = 0; i < fftObj->NN; i += 2)
     {
         fftObj->FFT_P[i + 1] = (fftObj->FFT_P[i + 0] <<= shift) + 1;
     }
 }

 /// <summary>
 /// Заполнение вектора изменения порядка следования данных перед FFT
 /// (для полифазного FFT)
 /// </summary>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void fill_FFT_PP(CFFT_Object *fftObj)
 {
     int i, j;
     
     // Выделяем память под вектор перестановок FFT...
     SAFE_DELETE(fftObj->FFT_PP);
     fftObj->FFT_PP = (int *)calloc(fftObj->NNPoly, sizeof(int));

     // Инициализация массива перестановок...
     memset(fftObj->FFT_PP, 0x00, (fftObj->NNPoly * sizeof(int)));

     // Заполняем вектор изменения порядка следования данных
     // (для полифазного FFT)...
     for(j = 0; j < LogX(fftObj->NPoly, 2); ++j)
     {
        for(i = 0; i < fftObj->NPoly; ++i)
        {
            fftObj->FFT_PP[i << 1] = ((fftObj->FFT_PP[i << 1] << 1) +
                                     ((i >> j) & 1));
        }
     }

     for(i = 0; i < fftObj->NNPoly; i += 2)
     {
         fftObj->FFT_PP[i + 1] = (fftObj->FFT_PP[i + 0] <<= 1) + 1;
     }
 }

 /// <summary>
 /// Заполнение вектора косинусного взвешивающего окна
 /// </summary>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void fill_FFT_TW_Cosine(CFFT_Object *fftObj)
 {     
     int i;

     // Формирующие параметры взвешивающего косинусного окна
     double arg, wval, a0, a1, a2, a3, ad;
     
     // Выделяем память под взвешивающее окно...
     fftObj->FFT_TW = (double *)calloc(fftObj->NN, sizeof(double));

     // Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)
     switch (fftObj->CosTW)
     {
        case RECTANGULAR_13dbPS:         { a0 = 1.0;       a1 = 0;         a2 = 0;         a3 = 0;         ad = 1.0;     break; }
        case HANN_31dbPS:                { a0 = 1.0;       a1 = 1.0;       a2 = 0;         a3 = 0;         ad = 2;       break; }
        case HAMMING_43dbPS:             { a0 = 0.54;      a1 = 0.46;      a2 = 0;         a3 = 0;         ad = 1.0;     break; }
        case MAX_ROLLOFF_3_TERM_46dbPS:  { a0 = 0.375;     a1 = 0.5;       a2 = 0.125;     a3 = 0;         ad = 1.0;     break; }
        case BLACKMAN_58dbPS:            { a0 = 0.42;      a1 = 0.5;       a2 = 0.08;      a3 = 0;         ad = 1.0;     break; }
        case COMPROMISE_3_TERM_64dbPS:   { a0 = 0.40897;   a1 = 0.5;       a2 = 0.09103;   a3 = 0;         ad = 1.0;     break; }
        case EXACT_BLACKMAN_68dbPS:      { a0 = 7938.0;    a1 = 9240.0;    a2 = 1430.0;    a3 = 0;         ad = 18608.0; break; }
        case MIN_SIDELOBE_3_TERM_71dbPS: { a0 = 0.4243801; a1 = 0.4973406; a2 = 0.0782793; a3 = 0;         ad = 1.0;     break; }
        case MAX_ROLLOFF_4_TERM_60dbPS:  { a0 = 10.0;      a1 = 15.0;      a2 = 6.0;       a3 = 1;         ad = 32.0;    break; }
        case COMPROMISE1_4_TERM_82dbPS:  { a0 = 0.338946;  a1 = 0.481973;  a2 = 0.161054;  a3 = 0.018027;  ad = 1.0;     break; }
        case COMPROMISE2_4_TERM_93dbPS:  { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
        default:
        case BLACKMAN_HARRIS_92dbPS:     { a0 = 0.35875;   a1 = 0.48829;   a2 = 0.14128;   a3 = 0.01168;   ad = 1.0;     break; }
        case NUTTALL_93dbPS:             { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
        case BLACKMAN_NUTTALL_98dbPS:    { a0 = 0.3635819; a1 = 0.4891775; a2 = 0.1365995; a3 = 0.0106411; ad = 1.0;     break; }
        case ROSENFIELD:                 { a0 = 0.762;     a1 = 1.0;       a2 = 0.238;     a3 = 0;         ad = a0;      break; }
     }

     // Заполняем взвешивающее окно коэффициентами...
     for(i = 0; i < fftObj->N; ++i)
     {
         arg  = (2.0 * M_PI * i) / (double)fftObj->N;
         wval = (a0 - a1 * cos(arg) + a2 * cos(2 * arg) - a3 * cos(3 * arg)) / ad;
         fftObj->FFT_TW[(i << 1) + 1] = fftObj->FFT_TW[(i << 1) + 0] = wval;
     }
 }

 /// <summary>
 /// Модифицир. функция Бесселя нулевого порядка первого рода
 /// </summary>
 /// <param name="arg"> Аргумент функции. </param>
 /// <returns> Значение функции. </returns>
 double BesselI0(double arg)
 {
     double numerator, denominator, z, z1, z2, z3, z4, z5,
            z6, z7, z8, z9, z10, z11, z12, z13, z_1, z_2;

     if(arg == 0.0)
     {
         return 1.0;
     }
     
     z = arg * arg;

     z1  = z * 0.210580722890567e-22 + 0.380715242345326e-19;
     z2  = z * z1  + 0.479440257548300e-16;
     z3  = z * z2  + 0.435125971262668e-13;
     z4  = z * z3  + 0.300931127112960e-10;
     z5  = z * z4  + 0.160224679395361e-7;
     z6  = z * z5  + 0.654858370096785e-5;
     z7  = z * z6  + 0.202591084143397e-2;
     z8  = z * z7  + 0.463076284721000e0;
     z9  = z * z8  + 0.754337328948189e2;
     z10 = z * z9  + 0.830792541809429e4;
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
 void fill_FFT_TW_Kaiser(CFFT_Object *fftObj)
 {
     int i, j;
     double norm, arg, w;

     // Выделяем память под вектор перестановок FFT...
     SAFE_DELETE(fftObj->FFT_TW);
     fftObj->FFT_TW = (double *)calloc(fftObj->NN, sizeof(double));

     // Нормирующий коэффициент окна Кайзера
     norm = BesselI0(fftObj->Beta);

     // Заполняем взвешивающее окно...
     for(i = 1; i <= (fftObj->N >> 1); ++i)
     {
         // arg = Beta * sqrt(1-(((2*(i-1))/(N-1))-1)^2);
         arg = fftObj->Beta *
                              sqrt(
                                    1 - pow(
                                             (
                                                (double)((i - 1) << 1)
                                              /
                                                (double)(fftObj->N - 1)
                                             ) - 1
                                        , 2)
                                   );

         w = BesselI0(arg) / norm;

         j = i - 1; // Приводим индекс от базы "1" к базе "0"
         fftObj->FFT_TW[(j << 1) + 0] = w; // left re
         fftObj->FFT_TW[(j << 1) + 1] = w; // left im
         fftObj->FFT_TW[(fftObj->NN - 2) - (j << 1) + 0] = w; // right re
         fftObj->FFT_TW[(fftObj->NN - 2) - (j << 1) + 1] = w; // right im
     }
 }

 /// <summary>
 /// "Контролер" объекта FFT
 /// </summary>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 /// <returns> Булевский флаг корректности конфигурации FFT. </returns>
 int CFFT_Inspector(CFFT_Object *fftObj)
 {
     // Размер окна не может быть меньше предельно-допустимого!
     if((fftObj->N     < MIN_FRAME_WIDTH) ||
        (fftObj->NPoly < MIN_FRAME_WIDTH) ||
        (fftObj->Beta  > MAX_KAISER_BETA))
     {
         return FALSE;
     }

     return TRUE;
 }

 /// <summary>
 /// "Деструктор" объекта FFT
 /// </summary>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void CFFT_Destructor(CFFT_Object *fftObj)
 {
     SAFE_DELETE(fftObj->FFT_P);
     SAFE_DELETE(fftObj->FFT_PP);
     SAFE_DELETE(fftObj->FFT_TW);
     SAFE_DELETE(fftObj);
 }

 /// <summary>
 /// Фабрика объектов FFT
 /// </summary>
 /// <param name="frameWidth"> Размер кадра. </param>
 /// <param name="cosTW"> Тип косинусного взвешивающего окна (если нет - используется окно Кайзера). </param>
 /// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
 /// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
 CFFT_Object * CFFT_Init(int frameWidth, int cosTW, double beta, int polyDiv2)
 {
     // Объект-результат
     CFFT_Object *fftObj = (CFFT_Object *)calloc(1, sizeof(CFFT_Object));

     // Заполнение полей объекта
     fftObj->N = ToLowerPowerOf2(frameWidth);  // Размер кадра FFT
     fftObj->NN = fftObj->N << 1;              // Кол-во точек (re + im)
     fftObj->NPoly = fftObj->N >> polyDiv2;    // Размер полифазного кадра FFT
     fftObj->CosTW = cosTW;                    // Тип косинусного взвешивающего окна
     fftObj->NNPoly = fftObj->NPoly << 1;      // Кол-во точек полифазного FFT
     fftObj->Beta = beta;                      // Форм-ий коэфф. окна Кайзера
     fftObj->PolyDiv = 1 << polyDiv2;          // Полифазный делитель

     fill_FFT_P(fftObj);  // Вектор изменения порядка след. данных перед FFT
     fill_FFT_PP(fftObj); // Вектор изменения порядка... (для полифазного FFT)
     
     if(fftObj->CosTW == NONE) //...если не задано взвешивающее окно косинусного типа
     {
         fill_FFT_TW_Kaiser(fftObj); // Взвешивающее окно Кайзера

     } else
     {
         fill_FFT_TW_Cosine(fftObj); // Косинусное взвешивающее окно
     }

     // Обрабатываем ситуацию со сбросом дампа...

#ifdef DUMP_MODE

     mkdir(DUMP_NAME);

     DumpInt(fftObj->FFT_P,     fftObj->NN,     DUMP_NAME, "FFT_P.int32");
     DumpInt(fftObj->FFT_PP,    fftObj->NNPoly, DUMP_NAME, "FFT_PP.int32");
     DumpDouble(fftObj->FFT_TW, fftObj->NN,     DUMP_NAME, "FFT_TW.double");

#endif

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
 CFFT_Object * CFFT_Constructor_Cosine(int frameWidth, int cosTW, int polyDiv2)
 {
     // Возвращаем объект "FFT"
     return CFFT_Init(frameWidth, cosTW, MAX_KAISER_BETA, polyDiv2);
 }

 /// <summary>
 /// Фабрика объектов FFT
 /// </summary>
 /// <param name="frameWidth"> Размер кадра. </param>
 /// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
 /// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
 CFFT_Object * CFFT_Constructor_Kaiser(int frameWidth, double beta, int polyDiv2)
 {
     // Возвращаем объект "FFT"
     return CFFT_Init(frameWidth, NONE, beta, polyDiv2);
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
 /// <param name="direction"> Направление преобразования (TRUE - прямое).
 /// </param>
 /// <param name="usePolyphase"> Использовать полифазное FFT? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void CFFT_Process(double *FFT_S, int FFT_S_Offset, double *FFT_T,
                   bool useTaperWindow, bool recoverAfterTaperWindow,
                   bool useNorm, bool direction, bool usePolyphase,
                   CFFT_Object *fftObj)
 {
     int i, j, mmax, isteps, n, istep, ii, m, jj;
     double isign, theta, sin05Th, wpr, wpi, wr, wi, tempr, tempi, wtemp;
     
     // Использование взвешивающего окна допустимо
     // только при прямом преобразовании
     if(direction && useTaperWindow)
     {
         if(!usePolyphase)
         {
             // Обычное FFT
             // только тогда, когда оно активно
             for(i = 0; i < fftObj->NN; ++i)
             {
                 FFT_T[i] = fftObj->FFT_TW[fftObj->FFT_P[i]] *
                            FFT_S[fftObj->FFT_P[i] + FFT_S_Offset];
             }
         }
         else
         {
             // Полифазное FFT
             // предполагает применение только на прямом проходе
             for(i = 0; i < fftObj->NNPoly; ++i)
             {
                 FFT_T[i] = 0;

                 // Накапливаем сумму текущей точки (в соответствии
                 // с количеством сегментов)
                 for(j = 0; j < fftObj->PolyDiv; ++j)
                 {
                     FFT_T[i] += fftObj->FFT_TW[fftObj->FFT_PP[i] +
                                 (j * fftObj->NNPoly)] * FFT_S[fftObj->FFT_PP[i] +
                                 (j * fftObj->NNPoly) + FFT_S_Offset];
                 }
             }
         }
     }
     else
     {
         // Обратный проход или прямой...
         // но без взвешивающего окна
         for(i = 0; i < fftObj->NN; ++i)
         {
             FFT_T[i] = FFT_S[fftObj->FFT_P[i] + FFT_S_Offset];
         }
     }

     // Нормализация коэффициентов производится при её выборе и только на
     // прямом проходе алгоритма (или если ситуация 100% симметрична)
     if((!direction) && (!useNorm))
     {
         for(i = 0; i < fftObj->NNPoly; ++i)
         {
             FFT_T[i] /= fftObj->N;
         }
     }

     // FFT Routine
     isign  = direction ? -1 : 1;
     mmax   = 2;
     isteps = 1;
     n = usePolyphase ? fftObj->NNPoly : fftObj->NN;
     while(n > mmax)
     {
         isteps++;
         istep   = mmax << 1;
         theta   = isign * ((2 * M_PI) / mmax);
         sin05Th = sin(0.5 * theta);
         wpr = -(2.0 * (sin05Th * sin05Th));
         wpi = sin(theta);
         wr  = 1.0;
         wi  = 0.0;

         for(ii = 1; ii <= (mmax >> 1); ++ii)
         {
             m = (ii << 1) - 1;
             for(jj = 0; jj <= ((n - m) >> isteps); ++jj)
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
             wr = wr * wpr - wi    * wpi + wr;
             wi = wi * wpr + wtemp * wpi + wi;
         }
         mmax = istep;
     }

     // Нормализация коэффициентов производится при её выборе и только
     // на прямом проходе алгоритма (или если ситуация 100% симметрична)
     if(direction && useNorm)
     {
         for(i = 0; i < n; ++i)
         {
             FFT_T[i] /= fftObj->N;
         }
     }

     // Аннигилируем взвешивающее окно (если оно равно нулю в некоторых точках
     // - результат восстановления неизвестен и полагаем его равным нулю)
     if((!direction) && useTaperWindow && recoverAfterTaperWindow)
     {
         for(i = 0; i < fftObj->NN; ++i)
         {
             FFT_T[i] = ((fftObj->FFT_TW[i] == 0) ?
                        0 : (FFT_T[i] / fftObj->FFT_TW[i]));
         }
     }
 }

 /// <summary>
 /// Исследование "левого" и "правого" каналов: ("левый" -
 /// действительная часть исходных данных, "правый" - мнимая часть)
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
 void CFFT_Explore(double *FFT_T, double *MagL, double *MagR, double *ACH,
                   double *ArgL, double *ArgR, double *PhaseLR, bool usePolyphase,
                   CFFT_Object *fftObj)
 {
     int N, i;
     double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
            lx, ly, rx, ry, argL, argR;

     // Определяем размерность FFT
     N = usePolyphase ? fftObj->NPoly : fftObj->N;

     // Вычисляем значения искомых величин для начальной точки
     // ("нулевая" гармоника)
     magL = FFT_T[0];
     magR = FFT_T[1];
          
     if(MagL    != NULL) MagL[0]    = magL;
     if(MagR    != NULL) MagR[0]    = magR;
     if(ACH     != NULL) ACH[0]     = magR / ((magL == 0) ? FLOAT_MIN : magL);
     if(ArgL    != NULL) ArgL[0]    = M_PI;
     if(ArgR    != NULL) ArgR[0]    = M_PI;
     if(PhaseLR != NULL) PhaseLR[0] = 0;

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

         magL = sqrt((lx * lx) + (ly * ly)) * 0.5;
         magR = sqrt((rx * rx) + (ry * ry)) * 0.5;
         argL = Safe_atan2(ly, lx);
         argR = Safe_atan2(ry, rx);

         if(MagL    != NULL) MagL[i] = magL;
         if(MagR    != NULL) MagR[i] = magR;
         if(ACH     != NULL) ACH[i]  = magR / ((magL == 0) ? FLOAT_MIN : magL);
         if(ArgL    != NULL) ArgL[i] = argL;
         if(ArgR    != NULL) ArgR[i] = argR;
         if(PhaseLR != NULL) PhaseLR[i] = PhaseNorm(argR - argL);
     }
 }

 /// <summary>
 /// Исследование результатов комплексного FFT (идентично CFFT из MathCAD)
 /// </summary>
 /// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
 /// <param name="Mag"> Магнитуды. </param>
 /// <param name="Arg"> Аргументы. </param>
 /// <param name="usePolyphase"> Использовать полифазное FFT? </param
 /// <param name="isMirror"> Зеркальное отображение спектра? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 void CFFT_ComplexExplore(double *FFT_T, double *Mag, double *Arg,
                          bool usePolyphase, bool isMirror,
                          CFFT_Object *fftObj)
 {
     int N, N_2, i;
     double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
            lx, ly, rx, ry, argL, argR;

     // Определяем размерность FFT
     N   = usePolyphase ? fftObj->NPoly : fftObj->N;
     N_2 = N >> 1;

     // Вычисляем значения искомых величин для начальной точки
     // ("нулевая" гармоника)
     lx = FFT_T[0];
     ly = FFT_T[1];
     rx = FFT_T[N + 0];
     ry = FFT_T[N + 1];
     
     if(Mag != NULL)
     {
         magL = sqrt((lx * lx) + (ly * ly));
         magR = sqrt((rx * rx) + (ry * ry));

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

     if(Arg != NULL)
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

         if(Mag != NULL)
         {
             magL = sqrt((lx * lx) + (ly * ly));
             magR = sqrt((rx * rx) + (ry * ry));

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

         if(Arg != NULL)
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
 void dB_Scale(double *data, double zero_db_level,
               CFFT_Object *fftObj)
 {
     int i;
     for(i = 0; i < fftObj->N >> 1; ++i)
     {
         data[i] = 10 * log10(data[i] / zero_db_level); // log
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
 CFFT_SelfTestResult SelfTest_S(double *FFT_S, double ACH_Difference,
                                CFFT_Object *fftObj)
 {
     double *FFT_S_backward, *FFT_T, *MagL, *MagR, *ACH, *ArgL, *ArgR, *PhaseLR;
     int N2, FFT_S_Offset, i;
     bool useTaperWindow, recoverAfterTaperWindow, useNorm, direction,
          usePolyphase;
     double maxDiff, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
            lx, ly, rx, ry, lx_, ly_, rx_, ry_, currentDiff;
     LARGE_INTEGER startCounter, CFFT_Process_counter, CFFT_Explore_counter, timerFrequency;
     int N_iters = 10000;

     // Cтруктура "Результат тестирования внутренней точности
     // прямого-обратного комплексного преобразования Фурье"
     CFFT_SelfTestResult selfTestResult;

     // Массив исходных данных - для заполнения на обратном ходе FFT
     FFT_S_backward = (double *)calloc(fftObj->NN, sizeof(double));

     // Целевой массив
     FFT_T = (double *)calloc(fftObj->NN, sizeof(double));

     // (Количество точек FFT / 2) - количество гармоник вместе с нулевой
     N2 = fftObj->N >> 1;

     // Массивы результатов Фурье-анализа
     MagL    = (double *)calloc(N2, sizeof(double));
     MagR    = (double *)calloc(N2, sizeof(double));
     ACH     = (double *)calloc(N2, sizeof(double));
     ArgL    = (double *)calloc(N2, sizeof(double));
     ArgR    = (double *)calloc(N2, sizeof(double));
     PhaseLR = (double *)calloc(N2, sizeof(double));

     // Не используем взвешивающее окно, но работаем
     // с нормализацией - направление прямое
     useTaperWindow = FALSE;
     FFT_S_Offset   = 0;
     recoverAfterTaperWindow = FALSE;
     useNorm      = TRUE;
     direction    = TRUE;
     usePolyphase = FALSE;
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
         FFT_T_N_i_Re = FFT_T[((fftObj->N - i) << 1) + 0];
         FFT_T_N_i_Im = FFT_T[((fftObj->N - i) << 1) + 1];

         lx = FFT_T_i_Re   + FFT_T_N_i_Re;
         ly = FFT_T_i_Im   - FFT_T_N_i_Im;
         rx = FFT_T_i_Im   + FFT_T_N_i_Im;
         ry = FFT_T_N_i_Re - FFT_T_i_Re;

         lx_ = 2 * MagL[i] * cos(ArgL[i]);
         ly_ = 2 * MagL[i] * sin(ArgL[i]);
         rx_ = 2 * MagR[i] * cos(ArgR[i]);
         ry_ = 2 * MagR[i] * sin(ArgR[i]);

         currentDiff = fabs(lx - lx_);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

         currentDiff = fabs(ly - ly_);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

         currentDiff = fabs(rx - rx_);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

         currentDiff = fabs(ry - ry_);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
     }

     // Сохраняем максимальную невязку перехода из алгебраической формы
     // в показательную и обратно
     selfTestResult.MaxDiff_ALG_to_EXP_to_ALG = maxDiff;

     // Не используем взвешивающее окно, но работаем
     // с нормализацией - направление обратное
     useTaperWindow = FALSE;
     FFT_S_Offset   = 0;
     recoverAfterTaperWindow = FALSE;
     useNorm      = TRUE;
     direction    = FALSE;
     usePolyphase = FALSE;
     CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
                  recoverAfterTaperWindow, useNorm, direction,
                  usePolyphase, fftObj);

     maxDiff = 0;
     for(i = 0; i < fftObj->N; ++i)
     {
         currentDiff = fabs(FFT_S_backward[i] - FFT_S[i]);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
     }

     // Сохраняем максимальную невязку после прямого-обратного преобразования
     // (без взвешивающего окна)
     selfTestResult.MaxDiff_FORWARD_BACKWARD = maxDiff;

     // Используем взвешивающее окно, работаем
     // с нормализацией - направление прямое
     useTaperWindow = TRUE;
     FFT_S_Offset   = 0;
     recoverAfterTaperWindow = FALSE;
     useNorm      = TRUE;
     direction    = TRUE;
     usePolyphase = FALSE;
     CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                  recoverAfterTaperWindow, useNorm, direction,
                  usePolyphase, fftObj);

     // Используем аннигиляцию взвешивающего окна, работаем
     // с нормализацией - направление обратное
     useTaperWindow = TRUE;
     FFT_S_Offset   = 0;
     recoverAfterTaperWindow = TRUE;
     useNorm      = TRUE;
     direction    = FALSE;
     usePolyphase = FALSE;
     CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
                  recoverAfterTaperWindow, useNorm, direction,
                  usePolyphase, fftObj);

     maxDiff = 0;
     for(i = (fftObj->NN / 2); i <= ((fftObj->NN * 3) / 4); ++i)
     {
         currentDiff = fabs(FFT_S_backward[i] - FFT_S[i]);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
     }

     // Сохраняем максимальную невязку после прямого-обратного
     // преобразования (c аннигиляцией взвешивающего окна)
     selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW = maxDiff;

     maxDiff = 0;
     for(i = 0; i < N2; ++i)
     {
         currentDiff = fabs(ACH[i] - ACH_Difference);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
     }

     // Сохраняем максимальную невязку по расчету заданной АЧХ
     selfTestResult.MaxDiff_ACH = maxDiff;

     maxDiff = 0;
     for(i = 0; i < N2; ++i)
     {
         currentDiff = fabs(PhaseLR[i]);
         maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
     }

     // Сохраняем максимальную невязку по расчету разности хода фаз каналов
     selfTestResult.MaxDiff_PhaseLR = maxDiff;

     // Performance Test
     QueryPerformanceFrequency(&timerFrequency);
     CFFT_Process_counter.QuadPart = CFFT_Explore_counter.QuadPart = 0;
               
     // Не используем взвешивающее окно, но работаем
     // с нормализацией - направление прямое
     useTaperWindow = FALSE;
     FFT_S_Offset   = 0;
     recoverAfterTaperWindow = FALSE;
     useNorm      = FALSE;
     direction    = TRUE;
     usePolyphase = FALSE;

     // CFFT_Process_time
     startCounter.QuadPart = 0;
     QueryPerformanceCounter(&startCounter);
     for (i = 0; i < N_iters; ++i)
     {
        CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                     recoverAfterTaperWindow, useNorm, direction,
                     usePolyphase, fftObj);
     }
     QueryPerformanceCounter(&CFFT_Process_counter);
     CFFT_Process_counter.QuadPart -= startCounter.QuadPart;
     selfTestResult.CFFT_Process_time  = (long double)CFFT_Process_counter.QuadPart / (long double)timerFrequency.QuadPart;
     selfTestResult.CFFT_Process_time /= (double)N_iters;

     // CFFT_Explore_time
     startCounter.QuadPart = 0;
     QueryPerformanceCounter(&startCounter);
     for (i = 0; i < N_iters; ++i)
     {
        CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                     usePolyphase, fftObj);
     }
     QueryPerformanceCounter(&CFFT_Explore_counter);
     CFFT_Explore_counter.QuadPart -= startCounter.QuadPart;
     selfTestResult.CFFT_Explore_time  = (long double)CFFT_Explore_counter.QuadPart / (long double)timerFrequency.QuadPart;
     selfTestResult.CFFT_Explore_time /= (double)N_iters;
          
     // Высвобождаем ресурсы динамической памяти
     SAFE_DELETE(FFT_S);
     SAFE_DELETE(FFT_S_backward);
     SAFE_DELETE(FFT_T);
     SAFE_DELETE(MagL);
     SAFE_DELETE(MagR);
     SAFE_DELETE(ACH);
     SAFE_DELETE(ArgL);
     SAFE_DELETE(ArgR);
     SAFE_DELETE(PhaseLR);

     // Проверка на допустимость полученных погрешностей
     if(selfTestResult.MaxDiff_ACH                     <= MAX_FFT_DIFF &&
        selfTestResult.MaxDiff_ALG_to_EXP_to_ALG       <= MAX_FFT_DIFF &&
        selfTestResult.MaxDiff_FORWARD_BACKWARD        <= MAX_FFT_DIFF &&
        selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW <= MAX_FFT_DIFF &&
        selfTestResult.MaxDiff_PhaseLR                 <= MAX_FFT_DIFF)
     {
         selfTestResult.AllOK = TRUE;

     } else
     {
         selfTestResult.AllOK = FALSE;
     }

#ifdef DUMP_MODE

     // Вердикт по точности выполнения самодиагностики
     DumpInt(&selfTestResult.AllOK,
             1, DUMP_NAME, "AllOK.int32");

     // Максимальная невязка по расчету заданной АЧХ
     DumpDouble(&selfTestResult.MaxDiff_ACH,
                1, DUMP_NAME, "MaxDiff_ACH.double");

     // Max. невязка ALG -> EXP и обратно
     DumpDouble(&selfTestResult.MaxDiff_ALG_to_EXP_to_ALG,
                1, DUMP_NAME, "MaxDiff_ALG_to_EXP_to_ALG.double");

     // Max. невязка FORVARD + BACKWARD
     DumpDouble(&selfTestResult.MaxDiff_FORWARD_BACKWARD,
                1, DUMP_NAME, "MaxDiff_FORWARD_BACKWARD.double");

     //...то же + восст. после TW
     DumpDouble(&selfTestResult.MaxDiff_FORWARD_BACKWARD_AntiTW,
                1, DUMP_NAME, "MaxDiff_FORWARD_BACKWARD_AntiTW.double");

     // Макс. невязка по расчету разности хода фаз
     DumpDouble(&selfTestResult.MaxDiff_PhaseLR,
                1, DUMP_NAME, "MaxDiff_PhaseLR.double");

     // Время работы CFFT_Process()
     DumpDouble(&selfTestResult.CFFT_Process_time,
                1, DUMP_NAME, "CFFT_Process_time.double");
     
     // Время работы CFFT_Explore()
     DumpDouble(&selfTestResult.CFFT_Explore_time,
                1, DUMP_NAME, "CFFT_Explore_time.double");

#endif

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
 CFFT_SelfTestResult SelfTest_RND(double ACH_Difference,
                                  CFFT_Object *fftObj)
 {    
     double randMult, randomValue;
     int i;
     
     // Массив исходных данных - стартовый
     double *FFT_S = (double *)calloc(fftObj->NN, sizeof(double));

     // Инициализируем генератор случайных чисел
     srand(time(NULL));

     // Множитель случайной выборки
     randMult = 1E07;

     // Заполняем исходный массив случайными данными...
     for(i = 0; i < fftObj->N; ++i)
     {
         // Получаем значение с датчика случайных чисел
         //(эмулируем ввод белого шума)...
         randomValue = (((double)rand() / (double)RAND_MAX) * randMult) -
                       (((double)rand() / (double)RAND_MAX) * randMult);

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
 int GetMaxIdx(double *data, int startIdx, int finishIdx)
 {
     int i, maxValIdxL, maxValIdxR;
     double currVal, maxVal;

     // Стартуем с первоначальным предположением...
     maxValIdxL = maxValIdxR = startIdx;
     maxVal = data[maxValIdxL];

     for(i = startIdx + 1; i <= finishIdx; ++i)
     {
         currVal = data[i];

         // Если текущее больше максимума -
         // сдвигаем оба индекса вправо...
         if(currVal > maxVal)
         {
             maxValIdxL = maxValIdxR = i;
             maxVal = currVal;

         } else
         //...а при равенстве максимуму
         // сдвигается только правый индекс
         if(currVal == maxVal)
         {
             maxValIdxR = i;
         }
     }

     return (maxValIdxL + maxValIdxR) >> 1;
 }
 
 /// <summary>
 /// Метод точного вычисления частоты по множеству гармоник
 /// </summary>
 /// <param name="Mag">  Магнитуды. </param>
 /// <param name="L">  Левая включенная граница для анализа. </param>
 /// <param name="R">  Правая включенная граница для анализа. </param>
 /// <param name="sampFreq"> Частота семплирования. </param>
 /// <param name="isComplex"> Комплексный режим? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
 double CalcExactFreq(double *Mag, int L, int R,
                      double sampFreq, bool isComplex,
                      CFFT_Object *fftObj)
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
     return FreqNode(exactFreqIdx, sampFreq, isComplex, fftObj);
 }

 /// <summary>
 /// Метод точного вычисления частоты по множеству гармоник
 /// </summary>
 /// <param name="Mag">  Магнитуды. </param>
 /// <param name="L">  Левая включенная граница для анализа. </param>
 /// <param name="R">  Правая включенная граница для анализа. </param>
 /// <param name="depth"> Глубина поиска. </param>
 /// <param name="sampFreq"> Частота семплирования. </param>
 /// <param name="isComplex"> Комплексный режим? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
 double ExactFreq(double *Mag, int L, int R, int depth,
                  double sampFreq, bool isComplex,
                  CFFT_Object *fftObj)
 {
     int fft_N, FFT_NodeMax, startIdx, finishIdx, startIdx2, finishIdx2,
         deltaStart, deltaFinish, deltaCorr;

     // Ищем максимальный индекс в установленных границах
     FFT_NodeMax = GetMaxIdx(Mag, L, R);
   
     // Вычисляем первоначальные индексы
     startIdx  = L;
     finishIdx = R;

     // Количество гармоник FFT
     fft_N = fftObj->N;

     // Корректируем индексы
     startIdx2  = startIdx  >= 0 ? startIdx : 0;
     finishIdx2 = finishIdx <= (fft_N - 1) ? finishIdx : (fft_N - 1);

     // Глубину обработки корректируем на максимальную дельту индексов
     deltaStart  = abs(startIdx2  - startIdx);
     deltaFinish = abs(finishIdx2 - finishIdx);
     deltaCorr   = max(deltaStart, deltaFinish);
     depth      -= deltaCorr;

     // Вычисляем окончательные индексы
     startIdx  = FFT_NodeMax - depth;
     finishIdx = FFT_NodeMax + depth;

     // Возвращаем точную частоту
     return CalcExactFreq(Mag, startIdx, finishIdx, sampFreq, isComplex, fftObj);
 }

 /// <summary>
 /// Метод точного вычисления частоты по множеству гармоник
 /// </summary>
 /// <param name="Mag">  Магнитуды. </param>
 /// <param name="depth"> Глубина поиска. </param>
 /// <param name="sampFreq"> Частота семплирования. </param>
 /// <param name="isComplex"> Комплексный режим? </param>
 /// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
 /// <returns> Точная частота, вычисленная по множеству гармоник. </returns>
 double ExactFreqAuto(double *Mag, int depth,
                      double sampFreq, bool isComplex,
                      CFFT_Object *fftObj)
 {
     return ExactFreq(Mag, 1, (fftObj->N >> 1) - 1, depth, sampFreq, isComplex, fftObj);
 }

#endif