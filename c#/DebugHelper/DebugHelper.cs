/*----------------------------------------------------------------------+
 |  filename:   DebugHelper.cs                                          |
 |----------------------------------------------------------------------|
 |  version:    1.00                                                    |
 |  revision:   24.12.2014  13:48                                       |
 |  author:     Дробанов Артём Федорович (DrAF)                         |
 |  e-mail:     draf@mail.ru                                            |
 |  purpose:    Класс-утилита для организации отладочного вывода        |
 |----------------------------------------------------------------------*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace DrAF.Utilities
{
    /// <summary>
    /// Отладочный хелпер
    /// </summary>
    public static class DebugHelper
    {
        /// <summary>
        /// Считывание значения типа "int"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Значение типа "int". </returns>
        public static int ReadInt(string dataPath, string fileName)
        {
            return ReadInts(dataPath, fileName)[0];
        }

        /// <summary>
        /// Считывание массива типа "int"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Массив типа "int". </returns>
        public static int[] ReadInts(string dataPath, string fileName)
        {
            var br = new BinaryReader(File.Open(Path.Combine(dataPath, fileName), FileMode.Open));
            var outList = new List<int>();

            try
            {
                while(true)
                {
                    outList.Add(br.ReadInt32());
                }
            }
            catch(EndOfStreamException)
            {
                br.Close();
            }

            if(outList.Count == 0)
            {
                throw new Exception("DebugHelper::ReadInts(): (outList.Count == 0)");
            }

            return outList.ToArray();
        }

        /// <summary>
        /// Запись значения типа "int"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteInt(string dataPath, string fileName, int data)
        {
            return WriteInts(dataPath, fileName, new int[] { data });
        }

        /// <summary>
        /// Запись последовательности типа "int"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteInts(string dataPath, string fileName, IEnumerable<int> data)
        {
            BinaryWriter bw = null;
            try
            {
                bw = new BinaryWriter(File.Open(Path.Combine(dataPath, fileName), FileMode.Create));

                foreach(var d in data)
                {
                    bw.Write(d);
                }

                bw.Flush();
                bw.Close();
            }
            catch(IOException)
            {
                if(bw != null)
                {
                    bw.Close();
                }
                return false;
            }

            return true;
        }

        /// <summary>
        /// Считывание значения типа "float"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Значение типа "float". </returns>
        public static float ReadFloat(string dataPath, string fileName)
        {
            return ReadFloats(dataPath, fileName)[0];
        }

        /// <summary>
        /// Считывание массива типа "float"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Массив типа "float". </returns>
        public static float[] ReadFloats(string dataPath, string fileName)
        {
            var br = new BinaryReader(File.Open(Path.Combine(dataPath, fileName), FileMode.Open));
            var outList = new List<float>();

            try
            {
                while (true)
                {
                    outList.Add(br.ReadSingle());
                }
            }
            catch(EndOfStreamException)
            {
                br.Close();
            }

            if(outList.Count == 0)
            {
                throw new Exception("DebugHelper::ReadFloats(): (outList.Count == 0)");
            }

            return outList.ToArray();
        }

        /// <summary>
        /// Запись значения типа "float"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteFloat(string dataPath, string fileName, float data)
        {
            return WriteFloats(dataPath, fileName, new float[] { data });
        }

        /// <summary>
        /// Запись последовательности типа "float"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteFloats(string dataPath, string fileName, IEnumerable<float> data)
        {
            BinaryWriter bw = null;
            try
            {
                bw = new BinaryWriter(File.Open(Path.Combine(dataPath, fileName), FileMode.Create));

                foreach(var d in data)
                {
                    bw.Write(d);
                }

                bw.Flush();
                bw.Close();
            }
            catch(IOException)
            {
                if(bw != null)
                {
                    bw.Close();
                }
                return false;
            }

            return true;
        }

        /// <summary>
        /// Считывание значения типа "double"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Значение типа "double". </returns>
        public static double ReadDouble(string dataPath, string fileName)
        {
            return ReadDoubles(dataPath, fileName)[0];
        }

        /// <summary>
        /// Считывание массива типа "double"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <returns> Массив типа "double". </returns>
        public static double[] ReadDoubles(string dataPath, string fileName)
        {
            var br = new BinaryReader(File.Open(Path.Combine(dataPath, fileName), FileMode.Open));
            var outList = new List<double>();

            try
            {
                while(true)
                {
                    outList.Add(br.ReadDouble());
                }
            }
            catch(EndOfStreamException)
            {
                br.Close();
            }

            if(outList.Count == 0)
            {
                throw new Exception("DebugHelper::ReadDoubles(): (outList.Count == 0)");
            }

            return outList.ToArray();
        }

        /// <summary>
        /// Запись значения типа "double"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteDouble(string dataPath, string fileName, double data)
        {
            return WriteDoubles(dataPath, fileName, new double[] { data });
        }

        /// <summary>
        /// Запись последовательности типа "double"
        /// </summary>
        /// <param name="dataPath"> Путь к данным. </param>
        /// <param name="fileName"> Имя файла дампа. </param>
        /// <param name="data"> Данные для сброса в дамп. </param>
        /// <returns> Булевский флаг операции. </returns>
        public static bool WriteDoubles(string dataPath, string fileName, IEnumerable<double> data)
        {
            BinaryWriter bw = null;
            try
            {
                bw = new BinaryWriter(File.Open(Path.Combine(dataPath, fileName), FileMode.Create));

                foreach(var d in data)
                {
                    bw.Write(d);
                }

                bw.Flush();
                bw.Close();
            }
            catch(IOException)
            {
                if(bw != null)
                {
                    bw.Close();
                }
                return false;
            }

            return true;
        }

        /// <summary>
        /// Метод "склеивания" двумерного массива в одномерный
        /// </summary>
        /// <typeparam name="T"> Тип данных массива. </typeparam>
        /// <param name="source"> Исходный массив. </param>
        /// <returns> Целевой массив. </returns>
        public static T[] Matrix2Vector<T>(T[][] source)
        {
            if(source == null)
            {
                throw new Exception("DebugHelper::Matrix2Vector<T>(): (source == null)");
            }

            if(source.Length == 0)
            {
                return new T[0];
            }

            if(source[0] == null)
            {
                throw new Exception("DebugHelper::Matrix2Vector<T>(): (source[0] == null)");
            }

            List<T> target = new List<T>(source.Length * source[0].Length);
            for(int i = 0; i < source.Length; i++)
            {
                // Извлекаем строку из источника и помещаем в приемник...
                target.AddRange(source[i]);

                // Подстраховка для доступа к массиву на следующей итерации...
                if((i < (source.Length - 1)) && (source[i + 1] == null))
                {
                    throw new Exception(String.Format("DebugHelper::Matrix2Vector<T>(): (source[{0}] == null)", i + 1));
                }
            }

            return target.ToArray();
        }
     
        /// <summary>
        /// Метод "склеивания" двумерного массива в одномерный
        /// </summary>
        /// <typeparam name="T"> Тип данных массива. </typeparam>
        /// <param name="source"> Исходный массив. </param>
        /// <returns> Целевой массив. </returns>
        public static T[] Matrix2Vector<T>(T[,] source)
        {
            if(source == null)
            {
                throw new Exception("DebugHelper::Matrix2Vector<T>(): (source == null)");
            }

            List<T> target = new List<T>(source.Length);
            for(int i = 0; i <= source.GetUpperBound(0); i++)
            {
                for(int j = 0; j <= source.GetUpperBound(1); j++)
                {
                    target.Add(source[i, j]);
                }
            }

            return target.ToArray();
        }

        /// <summary>
        /// Очистка массива
        /// </summary>
        /// <typeparam name="T"> Тип элементов массивов. </typeparam>
        /// <param name="array"> Массив для очистки. </param>
        public static void ClearArray<T>(T[] array)
        {
            if(array == null) return;
            Array.Clear(array, 0, array.Length);
        }

        /// <summary>
        /// Объединение двух массивов в результирующий
        /// </summary>
        /// <typeparam name="T"> Тип элементов массивов. </typeparam>
        /// <param name="array1"> Массив №1. </param>
        /// <param name="array2"> Массив №2. </param>
        /// <returns> Результирующий массив. </returns>
        public static T[] MergeArrays<T>(T[] array1, T[] array2)
        {
            if(array1 == null)
            {
                return (T[])array2.Clone();
            }

            if(array2 == null)
            {
                return (T[])array1.Clone();
            }

            var result = new T[array1.Length + array2.Length];

            Array.Copy(array1, 0, result, 0, array1.Length);
            Array.Copy(array2, 0, result, array1.Length, array2.Length);

            return result;
        }        
    }
}