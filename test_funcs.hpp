//
//  test_funcs.hpp
//  tacticsclone
//
//  Created by asdfuiop on 7/28/18.
//  Copyright © 2018 asdfuiop. All rights reserved.
//

#ifndef test_funcs_hpp
#define test_funcs_hpp

enum terrain {water, river, desert, scrub, tundra, boreal, temperate, wetland, tropical, tropical_dry, desert_high, temperate_high,
    boreal_high, tundra_high, wet_high, tropical_high, tropical_dry_high, dry_mount, snow_mount};

void checkthis(std::vector<std::vector<float>>& in, int rows, int cols, float min, float max);

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, double min, double max);

void checkthis(std::vector<std::vector<double>>& in, int rows, int cols, std::string filename, double min, double max);

void checkthis(std::vector<std::vector<terrain>>& in, int rows, int cols, std::string filename);

#endif /* test_funcs_hpp */
