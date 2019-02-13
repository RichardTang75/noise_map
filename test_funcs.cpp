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
enum terrain {water, river, desert, scrub, tundra, boreal, temperate, wetland, tropical, tropical_dry, desert_high, temperate_high,
    boreal_high, tundra_high, wet_high, tropical_high, tropical_dry_high, dry_mount, snow_mount};


void terrain_img(std::vector<std::vector<terrain>>& array, std::vector<unsigned char>& img,
                 const int& rows, const int& cols)
{
    std::unordered_map<terrain, std::vector<int>> terrain_colors;
    terrain_colors[terrain::water] = {92, 131, 249};
    terrain_colors[terrain::river] = {114, 184, 246};
    terrain_colors[terrain::tropical] = {175, 203, 91};
    terrain_colors[terrain::tropical_dry] = {203, 188, 92};
    terrain_colors[terrain::tropical_high] = {216, 240 ,142};
    terrain_colors[terrain::tropical_dry_high] = {215, 206, 134};
    terrain_colors[terrain::desert] = {236, 206, 98};
    terrain_colors[terrain::desert_high] = {255, 230, 139};
    terrain_colors[terrain::scrub] = {171, 225, 163};
    terrain_colors[terrain::boreal] = {50, 125, 107};
    terrain_colors[terrain::boreal_high] = {79, 154, 123};
    terrain_colors[terrain::tundra] = {93, 182, 144};
    terrain_colors[terrain::tundra_high] = {129, 209, 175};
    terrain_colors[terrain::temperate] = {84, 162, 74};
    terrain_colors[terrain::temperate_high] = {136, 204, 127};
    terrain_colors[terrain::wetland] = {65, 100, 61};
    terrain_colors[terrain::wet_high] = {99, 140, 94};
    terrain_colors[terrain::dry_mount] = {138, 128, 92};
    terrain_colors[terrain::snow_mount] = {249, 247, 237};
    for (int col = 0; col < cols; ++col)
    {
        for (int row = 0; row < rows; ++row)
        {
            int red = terrain_colors[array[col][row]][0];
            int green = terrain_colors[array[col][row]][1];
            int blue = terrain_colors[array[col][row]][2];
            img[4 * row*cols + 4 * col] = red;
            img[4 * row*cols + 4 * col + 1] = green;
            img[4 * row*cols + 4 * col + 2] = blue;
            img[4 * row*cols + 4 * col + 3] = 255;
        }
    }
}

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

void checkthis(std::vector<std::vector<terrain>>& in, int rows, int cols, std::string filename)
{
    std::vector<unsigned char> img(4 * rows*cols);
    terrain_img(in, img, rows, cols);
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
