/*----------------------------------------------------------------------+
 |  filename:   ExactFFT.h                                              |
 |----------------------------------------------------------------------|
 |  version:    7.10                                                    |
 |  revision:   07/09/2014  11:41                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Комплексное FFT                                         |
 |----------------------------------------------------------------------*/

 #ifndef _exactfft_h
 #define _exactfft_h

 #include <stdlib.h>
 #include <math.h>

 //------------------------------------
 //- ОТЛАДКА & DUMP
 //------------------------------------
 // !Закомментировать для отмены дампа
 #define DUMP_MODE
 #define DUMP_NAME "VS_C.dump"
 //------------------------------------

#ifdef DUMP_MODE

 #include <stdio.h>
 #include <direct.h>

#endif

 //------------------------
 //- БЛОК МАКРООПРЕДЕЛЕНИЙ
 //------------------------

 // Константы
 //                              3.14159265358979323846264338328
 #define M_PI                    3.14159265358979323846 //..264338328
 #define M_2PI                   2 * M_PI
 #define FLOAT_MIN               3.4E-38 // Допустимый минимум для операций float
 #define MAX_FFT_DIFF            1E-7    // Максимальная погрешность FFT
 #define MIN_FRAME_WIDTH         8       // Наименьший "рабочий" размер окна FFT
 #define MAX_KAISER_BETA         28      // Max. beta (Kaiser window), SL: ~240 dB
 #define MAX_PATH                256     // Максимальная длина пути

 // "Булевские" константы
 #define TRUE                    1       // "ИСТИНА"
 #define FALSE                   0       // "ЛОЖЬ"
 #define DIRECT                  1       // Обозначение прямого прохода FFT
 #define REVERSE                 0       // Обозначение обратного прохода FFT
 #define USING_NORM              1       // Используется нормализация
 #define NOT_USING_NORM          0       // Не используется нормализация
 #define USING_TAPER_WINDOW      1       // Используется взвешивающее окно
 #define NOT_USING_TAPER_WINDOW  0       // Не используется взвешивающее окно
 #define USING_POLYPHASE         1       // Используется полифазное FFT
 #define NOT_USING_POLYPHASE     0       // Не используется полифазное FFT

 // Безопасное высвобождение ресурсов указателя
 #define SAFE_DELETE(ptr)  if (ptr != NULL) \
                           { \
                               free(ptr); \
                               ptr = NULL; \
                           }
 
 // min / max
 #define min(a, b) ((a) < (b) ? (a) : (b))
 #define max(a, b) ((a) > (b) ? (a) : (b))

 // Вводим "булевский" тип данных
 typedef int bool;

 //---------------------------------------------
 //- Типы взвешивающих окон FFT
 //---------------------------------------------
 // Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)
 enum CosTW {
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

 //-------------------------
 //- Структура "Объект FFT"
 //-------------------------
 typedef struct
 {
     //-------------------------------------------------------------------------
     int     N;       // Количество точек FFT
     int     NN;      // Кол-во чисел (re + im) FFT
     int     NPoly;   // Пересчитанное для полифазного FFT количество точек
     int     NNPoly;  // Кол-во чисел (re + im) полифазного FFT
     int     CosTW;   // Тип косинусного взвешивающего окна (если нет - используется окно Кайзера)
     double  Beta;    // Формирующий коэффициент "beta" окна Кайзера     
     int     PolyDiv; // Делитель "полифазности" FFT ("0" - обычное FFT)
     //-------------------------------------------------------------------------
     int    *FFT_P;   // Вектор изменения порядка следования данных перед FFT
     int    *FFT_PP;  // Вектор изменения порядка... (для полифазного FFT)
     double *FFT_TW;  // Взвешивающее окно
     //-------------------------------------------------------------------------

 } CFFT_Object;

 //---------------------------------------------
 //- Структура "Результат самодиагностики CFFT"
 //---------------------------------------------
 typedef struct
 {
     //-------------------------------------------------------------------------
     bool AllOK; // Результат тестирования точности
     //-------------------------------------------------------------------------
     double MaxDiff_ACH; // Максимальная невязка по расчету заданной АЧХ
     double MaxDiff_ALG_to_EXP_to_ALG; // Max. невязка ALG -> EXP и обратно
     double MaxDiff_FORWARD_BACKWARD;  // Max. невязка FORVARD + BACKWARD
     double MaxDiff_FORWARD_BACKWARD_AntiTW; //...то же + восст. после TW
     double MaxDiff_PhaseLR;   // Макс. невязка по расчету разности хода фаз
     double CFFT_Process_time; // Время работы CFFT_Process()
     double CFFT_Explore_time; // Время работы CFFT_Explore()
     //-------------------------------------------------------------------------

 } CFFT_SelfTestResult;

#endif