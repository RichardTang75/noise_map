//
//  test_funcs.cpp
//  tacticsclone
//
//  Created by asdfuiop on 7/28/18.
//  Copyright © 2018 asdfuiop. All rights reserved.
//

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include "3rd party libs/lodepng/lodepng.h"


void array_img(std::vector<std::vector<double>>& array, std::vector<unsigned char>& img,
               const int& rows, const int& cols, double min, double max)
{
    //grayscale
    double range = max-min;
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            double normalized =(array[col][row]-min)/range;

            int gray = 255 * normalized;
            int blue =255*normalized;
            int red=255*normalized;
            int green=255*normalized;
//            if (normalized<.5)
//            {
//                red=0;
//                green=0;
//                blue=255-255*normalized;
//            }
//            else
//            {
//                red=255*normalized;
//                green=0;
//                blue=0;
//            }
//            if (normalized<.4)
//            {
//                red = 55+25*normalized;
//                green = 33+17*normalized;
//                blue = 230 + 25*normalized; //230-60
//            }
//            else if (normalized<.65)
//            {
//                blue = 25*normalized;
//                red=25*normalized;
//                green=170+25*normalized;
//            }
//            else if (normalized<.85)
//            {
//                red = 70*normalized;
//                blue = 10*normalized;
//                green = 90*normalized;
//            }
            img[4 * row*cols + 4 * col] = red;
            img[4 * row*cols + 4 * col + 1] = green;
            img[4 * row*cols + 4 * col + 2] = blue;
            img[4 * row*cols + 4 * col + 3] = 255;
        }
    }
}

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, double min, double max)
{
    ///set up the base
    //union it. set constraints on the other terrains.  (only add the initial if inside where you want)
    std::string filename = "ayoo.png";
    std::vector<unsigned char> img(4 * rows*cols);
    //vectormap map = set_to_vectormap(in, rows, cols);
    array_img(in, img, rows, cols, min, max);
    unsigned error = lodepng::encode(filename.c_str(), img, cols, rows);
    //    if (error)
    //    {
    //        std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << "\n";
    //    }
}

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, std::string filename, double min, double max)
{
    std::vector<unsigned char> img(4 * rows*cols);
    array_img(in, img, rows, cols, min, max);
    lodepng::encode(filename.c_str(), img, cols, rows);
}

void array_img(std::vector<std::vector<float>>& array, std::vector<unsigned char>& img,
               const int& rows, const int& cols, float min, float max)
{
    //grayscale
    float range = max-min;
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            float normalized =(array[col][row]-min)/range;
            int blue = 230-(normalized)*60; //70+35
            int red=55 + 25*normalized; //55+25
            int green=33 + 17*normalized; //33+17
            int gray = 241 - 241 * normalized;
            img[4 * row*cols + 4 * col] = red;
            img[4 * row*cols + 4 * col + 1] = green;
            img[4 * row*cols + 4 * col + 2] = blue;
            img[4 * row*cols + 4 * col + 3] = 255;
        }
    }
}

void checkthis(std::vector<std::vector<float>>& in, int rows, int cols, float min, float max)
{
    ///set up the base
    //union it. set constraints on the other terrains.  (only add the initial if inside where you want)
    std::string filename = "ayoo.png";
    std::vector<unsigned char> img(4 * rows*cols);
    //vectormap map = set_to_vectormap(in, rows, cols);
    array_img(in, img, rows, cols, min, max);
    unsigned error = lodepng::encode(filename.c_str(), img, rows, cols);
    //    if (error)
    //    {
    //        std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << "\n";
    //    }
}
