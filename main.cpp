//
//  main.cpp
//  fun with noise
//
//  Created by asdfuiop on 1/15/19.
//  Copyright © 2019 asdfuiop. All rights reserved.
//

//idea for river generation comes from https://cartographersguild.com/showthread.php?t=26931

#include <iostream>
#include <vector>
#include "test_funcs.hpp"
#include "3rd party libs/FastNoise/FastNoise.h"
#include "3rd party libs/poisson-disk-generator/PoissonGenerator.h"
#include "3rd party libs/ThreadPool/ThreadPool.h"
#include "3rd party libs/json/json.hpp"
#include "3rd party libs/delaunator-cpp/include/delaunator.hpp"
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <queue>
#include <fstream>
#include "river_creation.hpp"

enum class direction {north, south, east, west};

struct city_candidate
{
    coordinate location;
    double score;
};

struct city_and_nearby
{
    coordinate city;
    std::vector<coordinate> nearby;
};

struct bounds
{
    coordinate start;
    double radius;
};

struct city_flood_fill
{
    coordinate city_start;
    int number;
    std::unordered_map<coordinate, int, coordinate::hash> flood_fill;
    double max_rad; //purely for the second pass.
};

struct flood_fill_unit
{
    coordinate location;
    double cost;
    int manhattan_dist;
};
struct compare_flood_fill_units
{
    bool operator() (const flood_fill_unit& lhs, flood_fill_unit& rhs)
    {
        return lhs.manhattan_dist>rhs.manhattan_dist;
    }
};

bool compare_city_candidates(city_candidate& i, city_candidate& j)
{
    return i.score > j.score;
}


double dist_squared(int x1, int y1, int x2, int y2)
{
    return pow((x2 - x1), 2) + pow((y2 - y1), 2);
}



void draw_line(std::vector<std::vector<double>>& in, int startx, int starty, int endx, int endy)
{
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            in[next_x][next_y]=.25;
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            in[next_x][next_y]=.25;
        }
    }
}
//based off of http://members.chello.at/easyfilter/bresenham.html
void draw_line(std::vector<std::vector<double>>& in, double min, double max, int startx, int starty, int endx, int endy)
{
    int dx = abs(endx-startx);
    int dy = -abs(endy-starty);
    int stepx, stepy;
    if (endx>startx)
    {
        stepx=1;
    }
    else
    {
        stepx=-1;
    }
    if (endy>starty)
    {
        stepy=1;
    }
    else
    {
        stepy=-1;
    }
    int err = dx+dy;
    int twice_err;
    int next_x = startx;
    int next_y = starty;
    double normalized;
    while (true)
    {
        normalized = (in[next_x][next_y] - min) / (max-min);
        if (normalized>.4)
        {
            in[next_x][next_y]=min;
        }
        if (next_x == endx && next_y == endy)
        {
            break;
        }
        twice_err = 2*err;
        if (twice_err >= dy)
        {
            err+=dy;
            next_x+=stepx;
        }
        if (twice_err<=dx)
        {
            err+=dx;
            next_y+=stepy;
        }
    }
}
double raycast_cost(std::vector<std::vector<terrain>>& in, std::unordered_map<terrain, double>& terrain_costs, int startx, int starty, int endx, int endy)
{
    double total_cost=0;
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            total_cost+=terrain_costs[in[next_x][next_y]];
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            total_cost+=terrain_costs[in[next_x][next_y]];
        }
    }
    return total_cost;
}

void raycast_all_values_in_line(const std::vector<std::vector<terrain>>& in, std::unordered_map<terrain, double>& terrain_costs,
                                const int startx, const int starty, const int endx, const int endy,
                                std::unordered_map<coordinate, float, coordinate::hash>& costs)
{
    double total_cost=0;
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            total_cost+=terrain_costs[in[next_x][next_y]];
            costs[coordinate{next_x, next_y}]=total_cost;
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            total_cost+=terrain_costs[in[next_x][next_y]];
            costs[coordinate{next_x, next_y}]=total_cost;
        }
    }
}


void draw_line(std::unordered_set<coordinate, coordinate::hash> in, int startx, int starty, int endx, int endy)
{
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            in.emplace(coordinate{next_x, next_y});
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            in.emplace(coordinate{next_x, next_y});
        }
    }
}


void get_noise_maps(std::vector<std::vector<double>>& ridge_vectormap, std::vector<std::vector<double>>& moisture_vectormap, std::vector<std::vector<double>>& temperature_vectormap)
{
    FastNoise ridges;
    ridges.SetSeed(90823);
    ridges.SetNoiseType(FastNoise::SimplexFractal);
    ridges.SetFractalOctaves(8);
    ridges.SetFractalLacunarity(2);
    ridges.SetFractalGain(.6);
    ridges.SetFractalType(FastNoise::RigidMulti);
    ridges.SetGradientPerturbAmp(30);
    ridges.SetFrequency(.003);
    FastNoise elevation;
    elevation.SetSeed(283);
    elevation.SetNoiseType(FastNoise::SimplexFractal);
    elevation.SetFrequency(.003);
    elevation.SetFractalOctaves(8);
    FastNoise moisture;
    moisture.SetSeed(8743);
    moisture.SetNoiseType(FastNoise::SimplexFractal);
    moisture.SetFrequency(.007);
    moisture.SetGradientPerturbAmp(45);
    FastNoise temperature;
    temperature.SetSeed(2309);
    temperature.SetNoiseType(FastNoise::SimplexFractal);
    temperature.SetFrequency(.005);
    temperature.SetGradientPerturbAmp(45);
    for (int i=0; i<ridge_vectormap.size(); ++i)
    {
        std::vector<double> current_column = ridge_vectormap[i];
        for (int j=0; j<current_column.size(); ++j)
        {
            float perturb_x = float(i);
            float perturb_y = float(j);
            ridges.GradientPerturbFractal(perturb_x, perturb_y);
            ridge_vectormap[i][j] = .25*ridges.GetNoise(perturb_x, perturb_y) + elevation.GetNoise(perturb_x, perturb_y);
            float temperature_x = float(i);
            float temperature_y = float(j);
            temperature.GradientPerturb(temperature_x, temperature_y);
            temperature_vectormap[i][j] = temperature.GetNoise(temperature_x, temperature_y);
            float moisture_x = float(i);
            float moisture_y = float(j);
            moisture.GradientPerturb(moisture_x, moisture_y);
            moisture_vectormap[i][j] = moisture.GetNoise(temperature_x, temperature_y);
        }
    }
}
void find_max_min_vectormap(std::vector<std::vector<double>>& vectormap, double& min, double& max)
{
    double absolute_min=vectormap[0][0];
    double absolute_max=vectormap[0][0];
    for (std::vector<double> column : vectormap)
    {
        double current_min = *std::min_element(std::begin(column), std::end(column));
        if (current_min<absolute_min)
        {
            absolute_min=current_min;
        }
        double current_max = *std::max_element(std::begin(column), std::end(column));
        if(current_max>absolute_max)
        {
            absolute_max=current_max;
        }
    }
    min = absolute_min;
    max = absolute_max;
}

void find_max_min_vectormap(std::vector<std::vector<double>>& vectormap, double& min, double& max, double to_ignore)
{
    double absolute_min=vectormap[0][0];
    double absolute_max=vectormap[0][0];
    for (std::vector<double> column : vectormap)
    {
        double current_min = *std::min_element(std::begin(column), std::end(column));
        if ((current_min<absolute_min && current_min != to_ignore) || (absolute_min==to_ignore && current_min != to_ignore))
            {
                absolute_min=current_min;
            }
            double current_max = *std::max_element(std::begin(column), std::end(column));
            if(current_max>absolute_max)
            {
                absolute_max=current_max;
            }
            }
            min = absolute_min;
            max = absolute_max;
            }

//find all points_and_edges, separate the ones that are plausible, iterate through points
void get_edges_to_draw(std::vector<coordinate_edge>& edges, std::unordered_set<coordinate_edge, coordinate_edge::hash>& edges_to_draw,
                       std::vector<std::vector<double>>& ridge_vectormap, double absolute_min, double absolute_range)
{
    std::unordered_set<coordinate, coordinate::hash> visited_points;
    std::unordered_map<coordinate, point_and_edges, coordinate::hash> points_n_connections;
    for (coordinate_edge edge :  edges)
    {
        if (visited_points.count(edge.p1)==0)
        {
            point_and_edges temp = find_edges_involving_point(edges, edge.p1, ridge_vectormap, absolute_min, absolute_range);
            points_n_connections[temp.point] = temp;
            visited_points.emplace(edge.p1);
        }
    }
    std::unordered_map<coordinate, point_and_edges, coordinate::hash> coast_points;
    std::unordered_map<coordinate, point_and_edges, coordinate::hash> land_points;
    for (auto point : points_n_connections)
    {
        if (point.second.coastal)
        {
            coast_points[point.first] = point.second;
        }
        else if ((ridge_vectormap[int(point.first.x)][int(point.first.y)]-absolute_min)/(absolute_range) > .4)
        {
            land_points[point.first] = point.second;
        }
    }
    get_to_draw(coast_points, land_points, edges, edges_to_draw);
    //get_to_draw(temp.edge_and_probabilities, edges, edges_to_draw);
}

city_and_nearby get_nearby_cities(coordinate& city_loc, std::vector<coordinate_edge>& connections)
{
    city_and_nearby to_return;
    to_return.city =coordinate{int(city_loc.x), int(city_loc.y)};
    std::vector<coordinate> nearbys;
    for (coordinate_edge& connection :  connections)
    {
        if (connection.p1 == city_loc)
        {
            nearbys.push_back(coordinate{int(connection.p2.x), int(connection.p2.y)});
        }
        else if (connection.p2 == city_loc)
        {
            nearbys.push_back(coordinate{int(connection.p1.x), int(connection.p1.y)});
        }
    }
    to_return.nearby=nearbys;
    return to_return;
}
double get_raycast_max_radius(coordinate& city, std::vector<coordinate_edge>& all_edges)
{
    city_and_nearby nearby = get_nearby_cities(city, all_edges);
    int city_x = int(city.x);
    int city_y = int(city.y);
    double max_dist_squared=0;
    for (coordinate nearby_coord : nearby.nearby)
    {
        if (dist_squared(city_x, city_y, nearby_coord.x, nearby_coord.y)>max_dist_squared)
        {
            max_dist_squared = dist_squared(city_x, city_y, nearby_coord.x, nearby_coord.y);
        }
    }
    return max_dist_squared;
}

std::unordered_set<coordinate, coordinate::hash> second_flood_fill_consolidation (city_flood_fill& original,
                                                std::vector<std::vector<double>>& political_map,
                                                std::vector<std::vector<terrain>>& terrain_map, int num_rows, int num_cols)
{
    std::priority_queue<flood_fill_unit, std::vector<flood_fill_unit>, compare_flood_fill_units> possible;
    std::vector<coordinate> directions =
    {
        coordinate{0, -1},
        coordinate{-1, 0},                      coordinate{1, 0},
        coordinate{0, 1}
    };
    flood_fill_unit temp{original.city_start, 0, 0};
    possible.push(temp);
    std::unordered_set<coordinate, coordinate::hash> visited;
    std::unordered_set<coordinate, coordinate::hash> in_queue;
    while (possible.size()>0)
    {
        flood_fill_unit current = possible.top();
        possible.pop();
        visited.emplace(current.location);
        for (coordinate dir : directions)
        {
            coordinate next = coordinate{current.location.x+dir.x, current.location.y+dir.y};
            if (dist_squared(next.x, next.y, original.city_start.x, original.city_start.y)< original.max_rad &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0 && political_map[next.x][next.y] == original.number)
            {
                if (terrain_map[next.x][next.y] != terrain::water)
                {
                    possible.push(flood_fill_unit{next, 0, current.manhattan_dist+1});
                    in_queue.emplace(next);
                }
            }
        }
    }
    return visited;
}

city_flood_fill third_flood_fill_gaps (city_flood_fill& original, std::vector<std::vector<double>>& political_map,
                                       std::vector<std::vector<terrain>>& terrain_map, int num_rows, int num_cols)
{
    std::unordered_set<coordinate, coordinate::hash> contiguous_area = second_flood_fill_consolidation(original, political_map, terrain_map, num_rows, num_cols);
    city_flood_fill to_return;
    to_return.city_start = original.city_start;
    to_return.number = original.number;
    to_return.max_rad = original.max_rad;
    std::unordered_map<coordinate, int, coordinate::hash> costs;
    std::priority_queue<flood_fill_unit, std::vector<flood_fill_unit>, compare_flood_fill_units> possible;
    std::vector<coordinate> directions =
    {
        coordinate{0, -1},
        coordinate{-1, 0},                      coordinate{1, 0},
        coordinate{0, 1}
    };
    flood_fill_unit temp{original.city_start, 0, 0};
    possible.push(temp);
    costs[temp.location] = temp.cost;
    std::unordered_set<coordinate, coordinate::hash> visited;
    std::unordered_set<coordinate, coordinate::hash> in_queue;
    while (possible.size()>0)
    {
        flood_fill_unit current = possible.top();
        possible.pop();
        to_return.flood_fill[current.location] = current.cost;
        visited.emplace(current.location);
        for (coordinate dir : directions)
        {
            coordinate next = coordinate{current.location.x+dir.x, current.location.y+dir.y};
            if (dist_squared(next.x, next.y, original.city_start.x, original.city_start.y)< original.max_rad &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0)
            {
                if (terrain_map[next.x][next.y] != terrain::water)
                {
                    if (contiguous_area.count(next)>0)
                    {
                        possible.push(flood_fill_unit{next, 0, current.manhattan_dist+1});
                    }
                    else
                    {
                        possible.push(flood_fill_unit{next, current.cost+1, current.manhattan_dist+1});
                    }
                    in_queue.emplace(next);
                }
            }
        }
    }
    return to_return;
}

//do it in a circle
bool is_valid(const coordinate checking, const coordinate initial, const double max_rad_squared, const std::vector<coordinate>& directions)
{
    if (dist_squared(checking.x, checking.y, initial.x, initial.y) < max_rad_squared)
    {
        for (const coordinate& direction: directions)
        {
            if (dist_squared(checking.x+direction.x, checking.y+direction.y, initial.x, initial.y)>max_rad_squared)
            {
                return true;
            }
        }
    }
    return false;
}
std::unordered_set<coordinate, coordinate::hash> get_outer_bounds(double max_rad_square, coordinate initial)
{
    std::unordered_set<coordinate, coordinate::hash> to_return;
    int start_x = initial.x + std::ceil(sqrt(max_rad_square));
    int start_y = initial.y;
    if (dist_squared(start_x, start_y, initial.x, initial.y)>max_rad_square)
    {
        start_x-=1;
    }
    coordinate start{start_x, start_y};
    std::vector<coordinate> directions =
    {
        coordinate{-1, -1}, coordinate{0, -1}, coordinate{1, -1},
        coordinate{-1, 0},                      coordinate{1, 0},
        coordinate{-1, 1}, coordinate{0, 1}, coordinate{1, 1}
    };
    std::vector<coordinate> possible;
    possible.push_back(start);
    while (possible.size()>0)
    {
        coordinate current = possible.back();
        possible.pop_back();
        for (coordinate& direction : directions)
        {
            coordinate checking {current.x+direction.x, current.y+direction.y};
            if (is_valid(checking, initial, max_rad_square, directions))
            {
                possible.push_back(checking);
                to_return.emplace(checking);
            }
        }
    }
    return to_return;
}

city_flood_fill flood_fill_this_city (coordinate location, double max_rad_square, int number, std::vector<std::vector<terrain>>& map, int num_rows, int num_cols)
{
    city_flood_fill to_return;
    to_return.city_start = location;
    to_return.number = number;
    to_return.max_rad = max_rad_square;
    std::unordered_map<terrain, double> map_terrain_to_cost;
    map_terrain_to_cost[terrain::water] = 5;
    map_terrain_to_cost[terrain::river] = 40;
    map_terrain_to_cost[terrain::tropical] = 1;
    map_terrain_to_cost[terrain::tropical_dry] = 2;
    map_terrain_to_cost[terrain::tropical_high] = 1.5;
    map_terrain_to_cost[terrain::tropical_dry_high] = 3;
    map_terrain_to_cost[terrain::desert] = 4;
    map_terrain_to_cost[terrain::desert_high] = 6;
    map_terrain_to_cost[terrain::scrub] = 2;
    map_terrain_to_cost[terrain::boreal] = 1.5;
    map_terrain_to_cost[terrain::boreal_high] = 2;
    map_terrain_to_cost[terrain::tundra] = 3;
    map_terrain_to_cost[terrain::tundra_high] = 4.5;
    map_terrain_to_cost[terrain::temperate] = 1;
    map_terrain_to_cost[terrain::temperate_high] = 1.5;
    map_terrain_to_cost[terrain::wetland] = 1.5;
    map_terrain_to_cost[terrain::wet_high] = 2;
    map_terrain_to_cost[terrain::dry_mount] = 50;
    map_terrain_to_cost[terrain::snow_mount] = 50;
    std::unordered_map<coordinate, int, coordinate::hash> costs;
    std::vector<double> column(num_rows, 0);
    std::vector<std::vector<double>> what(num_cols, column);
    std::priority_queue<flood_fill_unit, std::vector<flood_fill_unit>, compare_flood_fill_units> possible;
    std::vector<coordinate> directions =
    {
        coordinate{-1, -1}, coordinate{0, -1}, coordinate{1, -1},
        coordinate{-1, 0},                      coordinate{1, 0},
        coordinate{-1, 1}, coordinate{0, 1}, coordinate{1, 1}
    };
    flood_fill_unit temp{location, 0, 0};
    possible.push(temp);
    costs[temp.location] = temp.cost;
    std::unordered_set<coordinate, coordinate::hash> visited;
    std::unordered_set<coordinate, coordinate::hash> in_queue;
    while (possible.size()>0)
    {
        flood_fill_unit current = possible.top();
        possible.pop();
        to_return.flood_fill[current.location] = current.cost /*+ pow(dist_squared(location.x, location.y, current.location.x, current.location.y), .5)*/;
        what[current.location.x][current.location.y] = current.cost /*+ pow(dist_squared(location.x, location.y, current.location.x, current.location.y), .5)*/;
        visited.emplace(current.location);
        for (coordinate dir : directions)
        {
            coordinate next = coordinate{current.location.x+dir.x, current.location.y+dir.y};
            if (dist_squared(next.x, next.y, location.x, location.y)<max_rad_square &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0)
            {
                if (map[next.x][next.y] != terrain::water)
                {
                    possible.push(flood_fill_unit{next, raycast_cost(map, map_terrain_to_cost, next.x, next.y, location.x, location.y), current.manhattan_dist+1});
                    in_queue.emplace(next);
                }
            }
        }
        
    }
    double what_min;
    double what_max;
    find_max_min_vectormap(what, what_min, what_max);
    what[location.x][location.y] =what_max;
    std::string file_name = "what"+std::to_string(number)+".png";
    checkthis(what, num_rows, num_cols, file_name, what_min, what_max);
    return to_return;
}
std::vector<city_candidate> rate_candidates(std::vector<coordinate>& poisson_city_points, std::vector<std::vector<terrain>>& terrain_vectormap, int radius, int num_rows, int num_columns)
{
    std::unordered_map<terrain, double> terrain_score;
    terrain_score[terrain::water] = 10;
    terrain_score[terrain::river] = 40;
    terrain_score[terrain::tropical] = 3;
    terrain_score[terrain::tropical_dry] = 1;
    terrain_score[terrain::tropical_high] = 2.5;
    terrain_score[terrain::tropical_dry_high] = .75;
    terrain_score[terrain::desert] = -2;
    terrain_score[terrain::desert_high] = -2.5;
    terrain_score[terrain::scrub] = -1;
    terrain_score[terrain::boreal] = 3;
    terrain_score[terrain::boreal_high] = 2.5;
    terrain_score[terrain::tundra] = -1;
    terrain_score[terrain::tundra_high] = -1.5;
    terrain_score[terrain::temperate] = 5;
    terrain_score[terrain::temperate_high] = 4;
    terrain_score[terrain::wetland] = 2;
    terrain_score[terrain::wet_high] = 2;
    terrain_score[terrain::dry_mount] = 0;
    terrain_score[terrain::snow_mount] = 0;
    std::vector<city_candidate> candidates;
    for (coordinate coord: poisson_city_points)
    {
        int x_value = coord.x;
        int y_value = coord.y;
        if (terrain_vectormap[x_value][y_value] != terrain::water) //river cities r cool
        {
            double score = 0;
            for (int i=x_value-radius; i<x_value+radius; i++)
            {
                for (int j=y_value-radius; j<y_value+radius; j++)
                {
                    if (i<num_columns && i>0 && j>0 && j<num_rows)
                    {
                        score+=terrain_score[terrain_vectormap[i][j]];
                    }
                }
            }
            candidates.push_back(city_candidate{coord, score});
        }
    }
    return candidates;
}

void fill_vectormap_with_cheapest(std::vector<std::vector<double>>& provincial_vectormap,
                                  std::vector<city_flood_fill>& city_flood_fills_for_provinces)
{
    for (int i=0; i<provincial_vectormap.size(); i++)
    {
        for (int j=0; j<provincial_vectormap[0].size(); j++)
        {
            int lowest_cost = -1;
            int lowest_number=-1;
            for (city_flood_fill& fill_n_cost : city_flood_fills_for_provinces)
            {
                if (fill_n_cost.flood_fill.count(coordinate{i, j}) > 0 )
                {
                    if (lowest_cost==-1)
                    {
                        lowest_cost = fill_n_cost.flood_fill[coordinate{i,j}];
                        lowest_number = fill_n_cost.number;
                    }
                    else if (lowest_cost > fill_n_cost.flood_fill[coordinate{i,j}])
                    {
                        lowest_cost = fill_n_cost.flood_fill[coordinate{i,j}];
                        lowest_number = fill_n_cost.number;
                    }
                }
            }
            provincial_vectormap[i][j] = double(lowest_number);
        }
    }
}
void urquhart_graph(std::vector<coordinate_triangle>& city_triangles, std::unordered_set<coordinate_edge, coordinate_edge::hash>& city_connections_not_draw)
{
    for (coordinate_triangle tri : city_triangles)
    {
        double len_tri_1 = dist_squared(tri.e1.p1.x, tri.e1.p1.y, tri.e1.p2.x, tri.e1.p2.y);
        double len_tri_2 = dist_squared(tri.e2.p1.x, tri.e2.p1.y, tri.e2.p2.x, tri.e2.p2.y);
        double len_tri_3 = dist_squared(tri.e3.p1.x, tri.e3.p1.y, tri.e3.p2.x, tri.e3.p2.y);
        std::vector<double> lengths = {len_tri_1, len_tri_2, len_tri_3};
        double longest_length =*std::max_element(lengths.begin(), lengths.end());
        if (longest_length == len_tri_1)
        {
            if (len_tri_1 != len_tri_2 && len_tri_1 != len_tri_3)
            {
                city_connections_not_draw.emplace(tri.e1);
            }
        }
        else if (longest_length == len_tri_2)
        {
            if (len_tri_1 != len_tri_2 && len_tri_2 != len_tri_3)
            {
                city_connections_not_draw.emplace(tri.e2);
            }
        }
        else
        {
            if (len_tri_2 != len_tri_3 && len_tri_1 != len_tri_3)
            {
                city_connections_not_draw.emplace(tri.e3);
            }
        }
    }
}

void drop (double& normalized, double& current_carrying, double& carrying_capacity, double& pickup_perc,
           double& perc_drop, double& min_drop_perc, int& current_x, int& current_y, std::vector<std::vector<double>>& to_effect)
{
    if (normalized < .4)
    {
        if (abs(current_carrying) < abs(carrying_capacity))
        {
            current_carrying += pickup_perc*carrying_capacity; //make moisture here max? doesn't really matter, it's water anyways
            if (abs(current_carrying) > abs(carrying_capacity))
            {
                current_carrying = carrying_capacity;
            }
        }
        /*                                              _________________________                                                   */
        to_effect[current_x][current_y] = 0;
        /*                                              _________________________                                                   */
    }
    else
    {
        double drop_amount = 0;
        if (current_carrying>0)
        {
            if (abs(current_carrying*perc_drop) > abs(min_drop_perc*carrying_capacity))
            {
                drop_amount = current_carrying*perc_drop;
            }
            else if (abs(current_carrying) > abs(min_drop_perc*carrying_capacity))
            {
                drop_amount = min_drop_perc*carrying_capacity;
            }
            else
            {
                drop_amount = current_carrying;
            }
            //rain-shadow attempt 1
            if (normalized>.65)
            {
                double diff = pow((normalized - .65)/.35, 2);
                if (abs(diff*carrying_capacity)>abs(drop_amount))
                {
                    current_carrying -= diff*carrying_capacity;
                }
                else
                {
                    current_carrying -= drop_amount;
                }
            }
            else
            {
                current_carrying -= drop_amount;
            }
        }
        to_effect[current_x][current_y] += drop_amount;
    }
}

void wind_effects (double carrying_capacity, double perc_drop, double pickup_perc, double min_drop_perc, direction wind_direction,
                   std::vector<std::vector<double>>& elevations, std::vector<std::vector<double>>& to_effect,
                   double absolute_min, double absolute_range)
{
    int current_x, current_y;
    if (wind_direction == direction::north || wind_direction == direction::south)
    {
        for (int i=0; i<elevations.size(); i++)
        {
            double current_carrying = 0;
            current_x = i;
            for (int j=0; j<elevations[0].size(); j++)
            {
                if (wind_direction == direction::north)
                {
                    current_y = j;
                }
                else
                {
                    current_y = elevations[0].size()-1-j;
                }
                double normalized = (elevations[current_x][current_y]-absolute_min)/absolute_range;
                drop(normalized, current_carrying, carrying_capacity, pickup_perc, perc_drop, min_drop_perc, current_x, current_y, to_effect);
            }
        }
    }
    else
    {
        for (int j=0; j<elevations[0].size(); j++)
        {
            double current_carrying = 0;
            current_y = j;
            for (int i=0; i<elevations.size(); i++)
            {
                if (wind_direction == direction::east)
                {
                    current_x = elevations.size()-1-i;
                }
                else
                {
                    current_x = i;
                }
                double normalized = (elevations[current_x][current_y]-absolute_min)/absolute_range;
                drop(normalized, current_carrying, carrying_capacity, pickup_perc, perc_drop, min_drop_perc, current_x, current_y, to_effect);
            }
        }
    }
}
nlohmann::json convert_vectormap(std::vector<std::vector<double>>& in)
{
    nlohmann::json equiv_vectormap;
    for (std::vector<double>& column: in)
    {
        nlohmann::json json_col(column);
        equiv_vectormap.push_back(json_col);
    }
    return equiv_vectormap;
}

std::vector<std::vector<double>> convert_json(nlohmann::json& in)
{
    std::vector<std::vector<double>> to_return;
    to_return = in.get<std::vector<std::vector<double>>>();
    return to_return;
}

nlohmann::json convert_to_json(std::vector<std::vector<terrain>>& in)
{
    nlohmann::json equiv_vectormap;
    for (std::vector<terrain>& column: in)
    {
        nlohmann::json json_col(column);
        equiv_vectormap.push_back(json_col);
    }
    return equiv_vectormap;
}

nlohmann::json convert_to_json(std::vector<std::vector<int>>& in)
{
    nlohmann::json equiv_vectormap;
    for (std::vector<int>& column: in)
    {
        nlohmann::json json_col(column);
        equiv_vectormap.push_back(json_col);
    }
    return equiv_vectormap;
}

std::vector<std::vector<terrain>> combined_vectormap(std::vector<std::vector<double>>& elevation, double& elev_min, double& elev_range,
                                                  std::vector<std::vector<double>>& moisture, double& moist_min, double& moist_range,
                                                  std::vector<std::vector<double>>& temperature, double& temperature_min, double& temperature_range)
{
    double water_height = .4;
    double mount_height = .85;
    double highland_height = .65;
    double high_temp = .75;
    double low_temp = .25;
    double fairly_dry = .3;
    double really_dry = .15;
    double fairly_wet = .7;
    std::vector<terrain> column (elevation[0].size());
    std::vector<std::vector<terrain>> final_terrain (elevation.size(), column);
    double normalized_elev, normalized_moist, normalized_temp;
    for (int i=0; i<elevation.size(); i++)
    {
        for (int j=0; j<elevation[0].size(); j++)
        {
            normalized_elev = (elevation[i][j]-elev_min) / elev_range;
            if (normalized_elev > water_height)
            {
                normalized_moist = (moisture[i][j]-moist_min) / moist_range;
                normalized_temp = (temperature[i][j]-temperature_min) / temperature_range;
                if (normalized_elev > mount_height)
                {
                    if (normalized_moist > fairly_dry)
                    {
                        final_terrain[i][j] = terrain::snow_mount;
                    }
                    else
                    {
                        final_terrain[i][j] = terrain::dry_mount;
                    }
                }
                else if (normalized_elev > highland_height)
                {
                    //tropical + tropical dry highlands?
                    if (normalized_temp < low_temp)
                    {
                        if (normalized_moist > fairly_dry)
                        {
                            final_terrain[i][j] = terrain::boreal_high;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::tundra_high;
                        }
                    }
                    else if (normalized_temp > high_temp)
                    {
                        if (normalized_moist > fairly_dry)
                        {
                            final_terrain[i][j] = terrain::tropical_high;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::tropical_dry_high;
                        }
                    }
                    else
                    {
                        if (normalized_moist > fairly_wet)
                        {
                            final_terrain[i][j] = terrain::wet_high;
                        }
                        else if (normalized_moist < really_dry)
                        {
                            final_terrain[i][j] = terrain::desert_high;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::temperate_high;
                        }
                    }
                }
                else
                {
                    /*scrub,temperate, tropical */
                    if (normalized_temp < low_temp)
                    {
                        if (normalized_moist > fairly_dry)
                        {
                            final_terrain[i][j] = terrain::boreal;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::tundra;
                        }
                    }
                    else if (normalized_temp > high_temp)
                    {
                        if (normalized_moist > fairly_dry)
                        {
                            final_terrain[i][j] = terrain::tropical;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::tropical_dry;
                        }
                    }
                    else
                    {
                        if (normalized_moist > fairly_wet)
                        {
                            final_terrain[i][j] = terrain::wetland;
                        }
                        else if (normalized_moist < really_dry)
                        {
                            final_terrain[i][j] = terrain::desert;
                        }
                        else if (normalized_moist < fairly_dry)
                        {
                            final_terrain[i][j] = terrain::scrub;
                        }
                        else
                        {
                            final_terrain[i][j] = terrain::temperate;
                        }
                    }
                }
            }
            else
            {
                if (normalized_elev==0)
                {
                    final_terrain[i][j] = terrain::river;
                }
                else
                {
                    final_terrain[i][j] = terrain::water;
                }
            }
        }
    }
    return final_terrain;
}

void edge_elevation_gradient(std::vector<std::vector<double>>& ridge_vectormap, int edge_size, double absolute_range,
                             int num_columns, int num_rows)
{
    for (int i=0; i<ridge_vectormap.size(); i++)
    {
        double column_edge_decrement=0;
        if (i<edge_size)
        {
            column_edge_decrement=(edge_size-i)*.9*absolute_range/edge_size;
        }
        else if (i>num_columns-edge_size)
        {
            column_edge_decrement=(i-(num_columns-edge_size))*.9*absolute_range/edge_size;
        }
        for (int j=0; j<ridge_vectormap[i].size(); j++)
        {
            double row_edge_decrement=0;
            if (j<edge_size)
            {
                row_edge_decrement=(edge_size-j)*.9*absolute_range/edge_size;
            }
            else if (j>num_rows-edge_size)
            {
                row_edge_decrement=(j-(num_rows-edge_size))*.9*absolute_range/edge_size;
            }
            double chosen_decrement = std::max(column_edge_decrement, row_edge_decrement);
            ridge_vectormap[i][j]=ridge_vectormap[i][j] - chosen_decrement;
        }
    }
}

int main(int argc, const char * argv[])
{
    int hard_conc = 1.5 * std::thread::hardware_concurrency();
    if (hard_conc==0)
    {
        hard_conc=3;
    }
    ThreadPool pool(hard_conc);
    
    bool generate_new_map = true;
    std::ifstream injson("elevations1.json");
    nlohmann::json elevation_json;
    if (injson.good())
    {
        generate_new_map=false;
        injson >> elevation_json;
    }
    
    int num_rows = 512; //1536;
    int num_columns = 512; //2560;
    std::vector<double> column(num_rows);
    std::vector<std::vector<double>> ridge_vectormap(num_columns, column);
    std::vector<std::vector<double>> moisture_vectormap(num_columns, column);
    std::vector<std::vector<double>> temperature_vectormap(num_columns, column);
    std::vector<std::vector<double>> poisson_vectormap(num_columns, column);
    std::vector<std::vector<double>> provincial_vectormap(num_columns, column);
    get_noise_maps(ridge_vectormap, moisture_vectormap, temperature_vectormap);
    double absolute_min;
    double absolute_max;
    find_max_min_vectormap(ridge_vectormap, absolute_min, absolute_max);
    std::cout<<absolute_min<<","<<absolute_max<<"\n";
    double absolute_range = (absolute_max-absolute_min);
    //go 20 in each side? set corners at -absolute, set 10 x 10 as half-absolute
    //add x+y, 40 = .9*absolute_range/40
    int edge_size = 160;
    edge_elevation_gradient(ridge_vectormap, edge_size, absolute_range,num_columns, num_rows);
    //don't recheck min and max before drawing
    double temperature_min;
    double temperature_max;
    find_max_min_vectormap(temperature_vectormap, temperature_min, temperature_max);
    double temperature_range = .8*(temperature_max-temperature_min);
    PoissonGenerator::DefaultPRNG some_generator = PoissonGenerator::DefaultPRNG(58);
    double width_to_height_ratio = double(num_columns)/double(num_rows);
    if (generate_new_map)
    {
        int num_points_desired = num_rows/4 * num_columns / 4;
        std::vector<PoissonGenerator::sPoint> poisson_points = PoissonGenerator::GeneratePoissonPoints(num_points_desired, some_generator, width_to_height_ratio, 30, false);
        std::vector<double> poisson_delaunay_points;
        for (auto i=poisson_points.begin(); i!=poisson_points.end(); i++)
        {
            int x_value = i->x*num_columns;
            int y_value = i->y*num_rows;
            poisson_delaunay_points.push_back(x_value);
            poisson_delaunay_points.push_back(y_value);
            poisson_vectormap[x_value][y_value] = .5;
        }
        
        delaunator::Delaunator river_triangulate(poisson_delaunay_points);
        std::vector<coordinate_edge> edges;
        for(std::size_t i = 0; i < river_triangulate.triangles.size(); i+=3)
        {
            int x0, x1, x2, y0, y1, y2;
            x0=river_triangulate.coords[2 * river_triangulate.triangles[i]];
            y0=river_triangulate.coords[2 * river_triangulate.triangles[i] + 1];
            x1=river_triangulate.coords[2 * river_triangulate.triangles[i + 1]];
            y1=river_triangulate.coords[2 * river_triangulate.triangles[i + 1] + 1];
            x2=river_triangulate.coords[2 * river_triangulate.triangles[i + 2]];
            y2=river_triangulate.coords[2 * river_triangulate.triangles[i + 2] + 1];
            coordinate p1 {x0, y0};
            coordinate p2 {x1, y1};
            coordinate p3 {x2, y2};
            coordinate_edge e1 {p1, p2};
            coordinate_edge e2 {p1, p3};
            coordinate_edge e3 {p2, p3};
            edges.push_back(e1);
            edges.push_back(e2);
            edges.push_back(e3);
        }
        std::unordered_set<coordinate_edge, coordinate_edge::hash> edges_to_draw;
        std::cout << "delaunay done" << edges.size()<<"\n";
        get_edges_to_draw(edges, edges_to_draw, ridge_vectormap, absolute_min, absolute_range);
//        std::future<void> please_save_time;
//        please_save_time = pool.enqueue(get_edges_to_draw, std::ref(edges), std::ref(edges_to_draw), std::ref(ridge_vectormap), absolute_min, absolute_range);
//        please_save_time.get(); //i swear i will actually make this save time when I clean it all up
        for (coordinate_edge edge_to_draw : edges_to_draw)
        {
            draw_line(poisson_vectormap, int(edge_to_draw.p1.x), int(edge_to_draw.p1.y),
                      (edge_to_draw.p2.x), int(edge_to_draw.p2.y));
            draw_line(ridge_vectormap, absolute_min, absolute_max, int(edge_to_draw.p1.x), int(edge_to_draw.p1.y),
                      (edge_to_draw.p2.x), int(edge_to_draw.p2.y));
        }
    }
    if (!generate_new_map)
    {
        ridge_vectormap = convert_json(elevation_json);
    }
    for (int i=0; i<temperature_vectormap.size(); i++)
    {
        for (int j=0; j<temperature_vectormap[i].size(); j++)
        {
            temperature_vectormap[i][j]-=j*2*temperature_range/num_rows;
            double normalized_elevation = ridge_vectormap[i][j]-absolute_min/(absolute_max-absolute_min);
            if (normalized_elevation>.75)
            {
                temperature_vectormap[i][j]-=(normalized_elevation-.75)*1.5*temperature_range;
            }
            else if (normalized_elevation<.4)
            {
                temperature_vectormap[i][j]-=.125*temperature_range;
            }
        }
    }
    double moisture_min, moisture_max;
    find_max_min_vectormap(moisture_vectormap, moisture_min, moisture_max);
    double moisture_range = moisture_max-moisture_min;
    //deal with rainshadows, perhaps run one without river-replenishment, remove rivers, then run again to refill rivers
    //160, 5, 30, 5
    wind_effects(moisture_range*16, .02, .2, .005, direction::east, ridge_vectormap, moisture_vectormap, absolute_min, absolute_range);
    wind_effects(moisture_range*.5, .16, .2, .04, direction::west, ridge_vectormap, moisture_vectormap, absolute_min, absolute_range);
    wind_effects(moisture_range*3, .08, .2, .01, direction::north, ridge_vectormap, moisture_vectormap, absolute_min, absolute_range);
    wind_effects(moisture_range*.5, .16, .2, .04, direction::south, ridge_vectormap, moisture_vectormap, absolute_min, absolute_range);
    wind_effects(-temperature_range*5, .02, .1, .005, direction::east, ridge_vectormap, temperature_vectormap, absolute_min, absolute_range);
    wind_effects(-temperature_range*1, .1, .1, .04, direction::west, ridge_vectormap, temperature_vectormap, absolute_min, absolute_range);
    wind_effects(-temperature_range*2, .05, .1, .01, direction::north, ridge_vectormap, temperature_vectormap, absolute_min, absolute_range);
    wind_effects(-temperature_range*1, .1, .1, .04, direction::south, ridge_vectormap, temperature_vectormap, absolute_min, absolute_range);
    PoissonGenerator::DefaultPRNG city_generator = PoissonGenerator::DefaultPRNG(85);
    find_max_min_vectormap(temperature_vectormap, temperature_min, temperature_max, 0);
    find_max_min_vectormap(moisture_vectormap, moisture_min, moisture_max, 0);
    temperature_range = temperature_max-temperature_min;
    moisture_range = moisture_max-moisture_min;
    std::cout<<moisture_min << "," << moisture_max <<"\n"<<temperature_min << "," << temperature_max << "\n";
    checkthis(ridge_vectormap, num_rows, num_columns, "ayoo3.png", absolute_min, absolute_max);
    checkthis(moisture_vectormap, num_rows, num_columns, "ayoo4.png", moisture_min, moisture_max);
    checkthis(temperature_vectormap, num_rows, num_columns, "ayoo5.png", temperature_min, temperature_max);
//    checkthis(poisson_vectormap, num_rows, num_columns, "ayoo14.png", 0, 1);
    std::vector<std::vector<terrain>> final_terrain_map = combined_vectormap(ridge_vectormap, absolute_min, absolute_range, moisture_vectormap, moisture_min, moisture_range, temperature_vectormap, temperature_min, temperature_range);
    checkthis(final_terrain_map, num_rows, num_columns, "finalmap1.png");
//    find_max_min_vectormap(ridge_vectormap, absolute_min, absolute_max);
//    std::cout<<absolute_min<<","<<absolute_max<<"\n";
    if (generate_new_map)
    {
        nlohmann::json elevations = convert_vectormap(ridge_vectormap);
        std::ofstream o("elevations1.json");
        o << elevations << std::endl;
    }
    
    /*------------------------------------------------------------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------------------------------*/
    /*------------------------------------------------Province Generation Here------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------------------------------------------------------------*/
    
    int num_city_sampling = num_rows/16 * num_columns / 16;
    std::vector<PoissonGenerator::sPoint> city_points = PoissonGenerator::GeneratePoissonPoints(num_city_sampling, city_generator, width_to_height_ratio, 30, false);
    std::vector<coordinate> poisson_city_points;
    for (auto i=city_points.begin(); i!=city_points.end(); i++)
    {
        int x_value = i->x*num_columns;
        int y_value = i->y*num_rows;
        poisson_city_points.push_back(coordinate{x_value, y_value});
    }
    int radius = 25;
    std::vector<city_candidate> candidates = rate_candidates(poisson_city_points, final_terrain_map, radius, num_rows, num_columns);
    std::sort(candidates.begin(), candidates.end(), compare_city_candidates);
    std::vector<city_candidate> chosen_candidates;
    bool good_candidate;
    int candidate_index = 0;
    while (chosen_candidates.size() < 25)
    {
        good_candidate=true;
        for (city_candidate already_chosen : chosen_candidates)
        {
            if (dist_squared(already_chosen.location.x, already_chosen.location.y, candidates[candidate_index].location.x, candidates[candidate_index].location.y)<pow(std::min(num_rows/8, num_columns/8), 2))
            {
                good_candidate=false;
                break;
            }
        }
        if (good_candidate)
        {
            chosen_candidates.push_back(candidates[candidate_index]);
        }
        candidate_index+=1;
        if (candidate_index>=candidates.size())
        {
            break;
        }
    }
    std::vector<double> city_poisson_points;
    for (city_candidate chosen: chosen_candidates)
    {
        poisson_vectormap[chosen.location.x][chosen.location.y]=1;
        city_poisson_points.push_back(chosen.location.x);
        city_poisson_points.push_back(chosen.location.y);
    }
    
    delaunator::Delaunator city_triangulate(city_poisson_points);
    std::vector<coordinate_edge> city_edges;
    std::vector<coordinate_triangle> city_triangles;
    for(std::size_t i = 0; i < city_triangulate.triangles.size(); i+=3)
    {
        int x0, x1, x2, y0, y1, y2;
        x0=city_triangulate.coords[2 * city_triangulate.triangles[i]];
        y0=city_triangulate.coords[2 * city_triangulate.triangles[i] + 1];
        x1=city_triangulate.coords[2 * city_triangulate.triangles[i + 1]];
        y1=city_triangulate.coords[2 * city_triangulate.triangles[i + 1] + 1];
        x2=city_triangulate.coords[2 * city_triangulate.triangles[i + 2]];
        y2=city_triangulate.coords[2 * city_triangulate.triangles[i + 2] + 1];
        coordinate p1 {x0, y0};
        coordinate p2 {x1, y1};
        coordinate p3 {x2, y2};
        coordinate_edge e1 {p1, p2};
        coordinate_edge e2 {p1, p3};
        coordinate_edge e3 {p2, p3};
        city_edges.push_back(e1);
        city_edges.push_back(e2);
        city_edges.push_back(e3);
        coordinate_triangle tri = {e1, e2, e3};
        city_triangles.push_back(tri);
    }
    std::unordered_map<coordinate, double, coordinate::hash> city_radius_bounds;
    for (city_candidate& chosen: chosen_candidates)
    {
         city_radius_bounds[coordinate{int(chosen.location.x), int(chosen.location.y)}] = get_raycast_max_radius(chosen.location, city_edges);
    }
    std::unordered_set<coordinate_edge, coordinate_edge::hash> city_connections_not_draw;
    urquhart_graph(city_triangles, city_connections_not_draw);
    std::vector<city_flood_fill> city_flood_fills_for_provinces;
    int number = 0;
    std::vector<std::future<city_flood_fill>> city_flood_fill_futures;

    for (city_candidate& chosen: chosen_candidates)
    {
        double max_rad_squared = city_radius_bounds[chosen.location];
        city_flood_fill_futures.push_back(pool.enqueue(flood_fill_this_city, chosen.location, max_rad_squared, number, std::ref(final_terrain_map), num_rows, num_columns));
        number+=1;
    }
    for (auto& future : city_flood_fill_futures)
    {
        city_flood_fills_for_provinces.push_back(future.get());
    }
    std::cout << "donezo";
    fill_vectormap_with_cheapest(provincial_vectormap, city_flood_fills_for_provinces);
    std::vector<city_flood_fill> finished_flood_fills;
    city_flood_fill_futures.clear();
    for (city_flood_fill& second_pass : city_flood_fills_for_provinces)
    {
        city_flood_fill_futures.push_back(pool.enqueue(third_flood_fill_gaps, std::ref(second_pass), provincial_vectormap, final_terrain_map, num_rows, num_columns));
    }
    for (auto& future : city_flood_fill_futures)
    {
        finished_flood_fills.push_back(future.get());
    }
    fill_vectormap_with_cheapest(provincial_vectormap, finished_flood_fills);
//    for (coordinate_edge city_edge : city_edges)
//    {
//        if (city_connections_not_draw.count(city_edge) == 0)
//        {
//            draw_line(provincial_vectormap, 0, 1, city_edge.p1.x, city_edge.p1.y, city_edge.p2.x, city_edge.p2.y);
//        }
//    }
    
    double provincial_min, provincial_max;
    find_max_min_vectormap(provincial_vectormap, provincial_min, provincial_max);
    checkthis(provincial_vectormap, num_rows, num_columns, "ayoo9.png", provincial_min, provincial_max);
    if (generate_new_map)
    {
        nlohmann::json terrains = convert_to_json(final_terrain_map);
        std::ofstream out_terrain("terrains.json");
        out_terrain << terrains << std::endl;
        nlohmann::json provinces = convert_vectormap(provincial_vectormap);
        std::ofstream out_province("provinces.json");
        out_province << provinces << std::endl;
    }
    return 0;
}
