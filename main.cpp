//
//  main.cpp
//  fun with noise
//
//  Created by asdfuiop on 1/15/19.
//  Copyright © 2019 asdfuiop. All rights reserved.
//

#include <iostream>
#include <vector>
#include "test_funcs.hpp"
#include "3rd party libs/FastNoise/FastNoise.h"
#include "3rd party libs/poisson-disk-generator/PoissonGenerator.h"
#include "3rd party libs/ThreadPool/ThreadPool.h"
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <queue>
#include "river_creation.hpp"

struct coordinate
{
    int x;
    int y;
    bool operator==(const coordinate& other) const
    {
        return (x==other.x && y==other.y);
    }
    struct hash
    {
        size_t operator() (const coordinate& coord) const
        {
            return coord.x*53+coord.y*31;
        }
    };
};

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
        for (int i=1; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            in[next_x][next_y]=.25;
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=1; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            in[next_x][next_y]=.25;
        }
    }
}
void draw_line(std::vector<std::vector<double>>& in, double min, double max, int startx, int starty, int endx, int endy)
{
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    double normalized = 0;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            normalized = (in[next_x][next_y] - min) / (max-min);
            if (normalized>.4)
            {
                in[next_x][next_y]=min;
            }
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            normalized = (in[next_x][next_y] - min) / (max-min);
            if (normalized>.4)
            {
                in[next_x][next_y]=min;
            }
        }
    }
}
double raycast_cost(std::vector<std::vector<double>>& in, double min, double range, int startx, int starty, int endx, int endy)
{
    double total_cost=0;
    int dx = endx - startx;
    int dy = endy - starty;
    int step;
    int next_x;
    int next_y;
    double normalized;
    if (abs(dx)>abs(dy))
    {
        step =dx/abs(dx);
        for (int i=0; i<abs(dx); i++)
        {
            next_x=i*step+startx;
            next_y=starty+i*dy/abs(dx);
            normalized = (in[next_x][next_y]-min) / range;
            if (normalized == 0)
            {
                total_cost+=7;
            }
            else if (normalized < .4)
            {
                total_cost+=3;
            }
            else if (normalized > .85)
            {
                total_cost+=5;
            }
            else if (normalized > .65)
            {
                total_cost+=2;
            }
            else
            {
                total_cost+=1;
            }
        }
    }
    else
    {
        step = dy/abs(dy);
        for (int i=0; i<abs(dy); i++)
        {
            next_x=startx+i*dx/abs(dy);
            next_y=i*step+starty;
            normalized = (in[next_x][next_y]-min) / range;
            if (normalized == 0)
            {
                total_cost+=7;
            }
            else if (normalized < .4)
            {
                total_cost+=3;
            }
            else if (normalized > .85)
            {
                total_cost+=5;
            }
            else if (normalized > .65)
            {
                total_cost+=2;
            }
            else
            {
                total_cost+=1;
            }
        }
    }
    return total_cost;
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
    ridges.SetGradientPerturbAmp(15);
    ridges.SetFrequency(.003);
    FastNoise elevation;
    elevation.SetSeed(283);
    elevation.SetNoiseType(FastNoise::SimplexFractal);
    elevation.SetFrequency(.003);
    elevation.SetFractalOctaves(8);
    FastNoise moisture;
    moisture.SetSeed(8743);
    moisture.SetNoiseType(FastNoise::SimplexFractal);
    moisture.SetFrequency(.003);
    ridges.SetGradientPerturbAmp(30);
    FastNoise temperature;
    temperature.SetSeed(2309);
    temperature.SetNoiseType(FastNoise::SimplexFractal);
    temperature.SetFrequency(.003);
    ridges.SetGradientPerturbAmp(30);
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


void get_edges_to_draw(std::vector<Edge<double>>& edges, std::unordered_set<Edge<double>, edge_hasher>& edges_to_draw,
                       std::vector<std::vector<double>>& ridge_vectormap, double absolute_min, double absolute_max)
{
    std::unordered_set<Vector2<double>, delaunay_point_hasher> visited_points;
    std::unordered_map<coordinate, point_and_edges, coordinate::hash> points_n_connections;
    for (Edge<double> edge :  edges)
    {
        if (visited_points.count(edge.p1)==0)
        {
            point_and_edges temp;
            temp.point = edge.p1;
            double normalized_start_elev = (ridge_vectormap[int(temp.point.x)][int(temp.point.y)]-absolute_min)/(absolute_max-absolute_min);
            temp.edge_and_probabilities = find_edges_involving_point(edges, temp.point, ridge_vectormap, absolute_min, absolute_max);
            get_to_draw(temp.edge_and_probabilities, edges, edges_to_draw);
            visited_points.emplace(edge.p1);
        }
    }
}

city_and_nearby get_nearby_cities(Vector2<double>& city_loc, std::vector<Edge<double>>& connections)
{
    city_and_nearby to_return;
    to_return.city =coordinate{int(city_loc.x), int(city_loc.y)};
    std::vector<coordinate> nearbys;
    for (Edge<double>& connection :  connections)
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
double get_raycast_max_radius(Vector2<double>& city, std::vector<Edge<double>>& all_edges)
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

std::unordered_set<coordinate, coordinate::hash> second_flood_fill_consolidation (city_flood_fill& original, std::vector<std::vector<double>>& political_map,
                                                std::vector<std::vector<double>>& elevation_map, double elevation_min,
                                                 double elevation_range, int num_rows, int num_cols)
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
            double normalized = (elevation_map[next.x][next.y]-elevation_min)/elevation_range;
            if (dist_squared(next.x, next.y, original.city_start.x, original.city_start.y)< original.max_rad &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0 && political_map[next.x][next.y] == original.number)
            {
                if (normalized > .4 || normalized == 0)
                {
                    possible.push(flood_fill_unit{next, 0, current.manhattan_dist+1});
                    in_queue.emplace(next);
                }
            }
        }
    }
    return visited;
}

city_flood_fill third_flood_fill_gaps (city_flood_fill& original, std::unordered_set<coordinate, coordinate::hash>& contiguous_area,
                                       std::vector<std::vector<double>>& political_map, std::vector<std::vector<double>>& elevation_map,
                                       double elevation_min, double elevation_range, int num_rows, int num_cols)
{
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
            double normalized = (elevation_map[next.x][next.y]-elevation_min)/elevation_range;
            if (dist_squared(next.x, next.y, original.city_start.x, original.city_start.y)< original.max_rad &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0)
            {
                if (normalized > .4 || normalized == 0)
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

city_flood_fill flood_fill_this_city (coordinate& location, double& max_rad_square, int number, std::vector<std::vector<double>>& map, double absolute_min, double absolute_range, int num_rows, int num_cols)
{
    city_flood_fill to_return;
    to_return.city_start = location;
    to_return.number = number;
    to_return.max_rad = max_rad_square;
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
        to_return.flood_fill[current.location] = current.cost + pow(dist_squared(location.x, location.y, current.location.x, current.location.y), .5);
        what[current.location.x][current.location.y] = current.cost + pow(dist_squared(location.x, location.y, current.location.x, current.location.y), .5);
        visited.emplace(current.location);
        for (coordinate dir : directions)
        {
            coordinate next = coordinate{current.location.x+dir.x, current.location.y+dir.y};
            double normalized = (map[next.x][next.y]-absolute_min)/absolute_range;
            if (dist_squared(next.x, next.y, location.x, location.y)<max_rad_square &&
                next.x > 0 && next.x < num_cols && next.y >0 && next.y < num_rows &&
                visited.count(next) == 0 && in_queue.count(next) == 0)
            {
                if (normalized > .4 || normalized == 0)
                {
                    possible.push(flood_fill_unit{next, raycast_cost(map, absolute_min, absolute_range, next.x, next.y, location.x, location.y), current.manhattan_dist+1});
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
    checkthis(what, 512, 512, file_name, what_min, what_max);
    return to_return;
}
std::vector<city_candidate> rate_candidates(std::vector<coordinate>& poisson_city_points, std::vector<std::vector<double>>& ridge_vectormap, double absolute_min, double absolute_range, int radius, int num_rows, int num_columns, std::vector<std::vector<double>>& temperature_vectormap, double temperature_min, double temperature_range)
{
    std::vector<city_candidate> candidates;
    for (coordinate coord: poisson_city_points)
    {
        int x_value = coord.x;
        int y_value = coord.y;
        if ((ridge_vectormap[x_value][y_value]-absolute_min)/absolute_range > .4)
        {
            double score = 0;
            for (int i=x_value-radius; i<x_value+radius; i++)
            {
                for (int j=y_value-radius; j<y_value+radius; j++)
                {
                    if (i<num_columns && i>0 && j>0 && j<num_rows)
                    {
                        if (ridge_vectormap[i][j]==absolute_min)
                        {
                            score+=25;
                        }
                        else if ((ridge_vectormap[i][j]-absolute_min)/absolute_range < .4)
                        {
                            score+=5;
                        }
                        else if ((ridge_vectormap[i][j]-absolute_min)/absolute_range < .5)
                        {
                            score+=1;
                        }
                        else if ((ridge_vectormap[i][j]-absolute_min)/absolute_range < .65)
                        {
                            score+=.5;
                        }
                        else if ((ridge_vectormap[i][j]-absolute_min)/absolute_range > .85)
                        {
                            score-=.5;
                        }
                        if ((temperature_vectormap[i][j]-temperature_min)/(temperature_range) < .2)
                        {
                            score-=.15;
                        }
                        else if ((temperature_vectormap[i][j]-temperature_min)/(temperature_range) > .8)
                        {
                            score-=.15;
                        }
                        else if ((temperature_vectormap[i][j]-temperature_min)/(temperature_range) > .4 && (temperature_vectormap[i][j]-temperature_min)/(temperature_range) < .6)
                        {
                            score+=.15;
                        }
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
            int lowest_number=11;
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
void urquhart_graph(std::vector<Triangle<double>>& city_triangles, std::unordered_set<Edge<double>, edge_hasher>& city_connections_not_draw)
{
    for (Triangle<double> tri : city_triangles)
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
int main(int argc, const char * argv[])
{
    int hard_conc = 1.5 * std::thread::hardware_concurrency();
    if (hard_conc==0)
    {
        hard_conc=3;
    }
    ThreadPool pool(hard_conc);
    
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
    //don't recheck min and max before drawing
    double temperature_min;
    double temperature_max;
    find_max_min_vectormap(temperature_vectormap, temperature_min, temperature_max);
    double temperature_range = .8*(temperature_max-temperature_min);
    double oceanic_cooling_power = temperature_range/8; //do i really want to do this, probably not
    PoissonGenerator::DefaultPRNG some_generator = PoissonGenerator::DefaultPRNG(58);
    int num_points_desired = num_rows/16 * num_columns / 16;
    double width_to_height_ratio = double(num_columns)/double(num_rows);
    std::vector<PoissonGenerator::sPoint> poisson_points = PoissonGenerator::GeneratePoissonPoints(num_points_desired, some_generator, width_to_height_ratio, 30, false);
    std::vector<Vector2<double>> poisson_delaunay_points;
    for (auto i=poisson_points.begin(); i!=poisson_points.end(); i++)
    {
        int x_value = i->x*num_columns;
        int y_value = i->y*num_rows;
        poisson_delaunay_points.push_back(Vector2<double>(x_value, y_value));
        poisson_vectormap[x_value][y_value] = .5;
    }
    Delaunay<double> triangulation;
    std::vector<Triangle<double> > triangles = triangulation.triangulate(poisson_delaunay_points);
    std::cout << triangles.size() << " triangles generated\n";
    std::vector<Edge<double> > edges = triangulation.getEdges();
    std::unordered_set<Edge<double>, edge_hasher> edges_to_draw;
    std::future<void> please_save_time;
    please_save_time = pool.enqueue(get_edges_to_draw, std::ref(edges), std::ref(edges_to_draw), std::ref(ridge_vectormap), absolute_min, absolute_max);
    please_save_time.get(); //i swear i will actually make this save time when I finish the province stuff
    for (Edge<double> edge_to_draw : edges_to_draw)
    {
        draw_line(poisson_vectormap, int(edge_to_draw.p1.x), int(edge_to_draw.p1.y),
                                            (edge_to_draw.p2.x), int(edge_to_draw.p2.y));
        draw_line(ridge_vectormap, absolute_min, absolute_max, int(edge_to_draw.p1.x), int(edge_to_draw.p1.y),
                  (edge_to_draw.p2.x), int(edge_to_draw.p2.y));
    }
    for (int i=0; i<temperature_vectormap.size(); i++)
    {
        for (int j=0; j<temperature_vectormap[i].size(); j++)
        {
            temperature_vectormap[i][j]-=j*temperature_range/num_rows;
            double normalized_elevation = ridge_vectormap[i][j]-absolute_min/(absolute_max-absolute_min);
            if (normalized_elevation>.75)
            {
                temperature_vectormap[i][j]-=(normalized_elevation-.75)*1.5*temperature_range;
            }
        }
    }
    PoissonGenerator::DefaultPRNG city_generator = PoissonGenerator::DefaultPRNG(85);
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
    std::vector<city_candidate> candidates = rate_candidates(poisson_city_points, ridge_vectormap, absolute_min, absolute_range, radius, num_rows, num_columns, temperature_vectormap, temperature_min, temperature_range);
    std::sort(candidates.begin(), candidates.end(), compare_city_candidates);
    std::vector<city_candidate> chosen_candidates;
    bool good_candidate;
    int candidate_index = 0;
    while (chosen_candidates.size() < 10)
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
    std::vector<Vector2<double>> points_for_bounds;
    Delaunay<double> city_bound_triangulation;
    for (city_candidate chosen: chosen_candidates)
    {
        poisson_vectormap[chosen.location.x][chosen.location.y]=1;
        points_for_bounds.push_back(Vector2<double>(chosen.location.x, chosen.location.y));
    }
    std::vector<Triangle<double> > city_triangles = city_bound_triangulation.triangulate(points_for_bounds);
    std::vector<Edge<double> > city_edges = city_bound_triangulation.getEdges();
    std::unordered_map<coordinate, double, coordinate::hash> city_radius_bounds;
    for (Vector2<double> chosen : points_for_bounds)
    {
        city_radius_bounds[coordinate{int(chosen.x), int(chosen.y)}] = get_raycast_max_radius(chosen, city_edges);
    }
    std::unordered_set<Edge<double>, edge_hasher> city_connections_not_draw;
    urquhart_graph(city_triangles, city_connections_not_draw);
    for (Edge<double> city_edge : city_edges)
    {
        if (city_connections_not_draw.count(city_edge) == 0)
        {
            draw_line(poisson_vectormap, city_edge.p1.x, city_edge.p1.y, city_edge.p2.x, city_edge.p2.y);
        }
        draw_line(poisson_vectormap, city_edge.p1.x, city_edge.p1.y, city_edge.p2.x, city_edge.p2.y);
    }
    std::vector<city_flood_fill> city_flood_fills_for_provinces;
    int number = 0;
    for (city_candidate chosen: chosen_candidates)
    {
        double max_rad_squared = city_radius_bounds[chosen.location];
        city_flood_fills_for_provinces.push_back(flood_fill_this_city(chosen.location, max_rad_squared, number, ridge_vectormap, absolute_min, absolute_range, num_rows, num_columns));
        number+=1;
    }
    fill_vectormap_with_cheapest(provincial_vectormap, city_flood_fills_for_provinces);
    std::vector<city_flood_fill> finished_flood_fills;
    for (city_flood_fill& second_pass : city_flood_fills_for_provinces)
    {
        std::unordered_set<coordinate, coordinate::hash> contiguous = second_flood_fill_consolidation(second_pass, provincial_vectormap, ridge_vectormap, absolute_min, absolute_range, num_rows, num_columns);
        finished_flood_fills.push_back(third_flood_fill_gaps(second_pass, contiguous, provincial_vectormap, ridge_vectormap, absolute_min, absolute_range, num_rows, num_columns));
    }
    fill_vectormap_with_cheapest(provincial_vectormap, finished_flood_fills);
    double moisture_min, moisture_max;
    find_max_min_vectormap(moisture_vectormap, moisture_min, moisture_max);
    find_max_min_vectormap(temperature_vectormap, temperature_min, temperature_max);
//    checkthis(ridge_vectormap, num_rows, num_columns, "ayoo8.png", absolute_min, absolute_max);
//    checkthis(moisture_vectormap, num_rows, num_columns, "ayoo1.png", moisture_min, moisture_max);
//    checkthis(temperature_vectormap, num_rows, num_columns, "ayoo2.png", temperature_min, temperature_max);
    checkthis(poisson_vectormap, num_rows, num_columns, "ayoo7.png", 0, 1);
    double provincial_min, provincial_max;
    find_max_min_vectormap(provincial_vectormap, provincial_min, provincial_max);
    checkthis(provincial_vectormap, num_rows, num_columns, "ayoo9.png", provincial_min, provincial_max);
    find_max_min_vectormap(ridge_vectormap, absolute_min, absolute_max);
    std::cout<<absolute_min<<","<<absolute_max<<"\n";
    temperature_vectormap[5][12]=5;
    std::cout << "Hello, World!\n";
    return 0;
}

