//
//  test_funcs.hpp
//  tacticsclone
//
//  Created by asdfuiop on 7/28/18.
//  Copyright © 2018 asdfuiop. All rights reserved.
//

#ifndef test_funcs_hpp
#define test_funcs_hpp

void checkthis(std::vector<std::vector<float>>& in, int rows, int cols, float min, float max);

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, double min, double max);

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, std::string filename, double min, double max);

#endif /* test_funcs_hpp */
