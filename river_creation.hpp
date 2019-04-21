//
//  river_creation.hpp
//  fun with noise
//
//  Created by asdfuiop on 2/5/19.
//  Copyright © 2019 asdfuiop. All rights reserved.
//

#ifndef river_creation_hpp
#define river_creation_hpp

#include <random>

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

struct coordinate_edge
{
    coordinate p1;
    coordinate p2;
    bool operator==(const coordinate_edge& other) const
    {
        return (p1==other.p1 && p2==other.p2);
    }
    struct hash
    {
        std::size_t operator()(const coordinate_edge& e) const
        {
            return (e.p1.x*53 + e.p1.y*31 + e.p2.x*53 + e.p2.y*31);
        }
    };
};

struct coordinate_triangle
{
    coordinate_edge e1;
    coordinate_edge e2;
    coordinate_edge e3;
};

/*
 struct edge_hasher
 {
 std::size_t operator()(Edge<double> const& e) const
 {
 return (e.p1.x*53 + e.p1.y*31 + e.p2.x*53 + e.p2.y*31);
 }
 };
 */

//should redo this to ensure no doubles, use unordered set.
struct edge_information
{
    int edge_index;
    double edge_probability;
};

struct growing_river
{
    std::vector<coordinate> grow_points;
    std::vector<coordinate> old_grow_points; //used for culling
    std::unordered_set<coordinate, coordinate::hash> all_points;
};

bool compare_edge_probabilities(edge_information& i, edge_information& j)
{
    return i.edge_probability > j.edge_probability;
}

struct point_and_edges
{
    coordinate point;
    std::vector<edge_information> edge_and_probabilities;
    bool coastal;
};

std::vector<int> find_edges_involving_point(std::vector<coordinate_edge>& edges, coordinate& point)
{
    std::vector<int> to_return_indices;
    for (int i=0; i<edges.size(); i++)
    {
        coordinate_edge checkin = edges[i];
        if ((checkin.p1 == point) || (checkin.p2 == point))
        {
            to_return_indices.push_back(i);
        }
    }
    return to_return_indices;
}

point_and_edges find_edges_involving_point(std::vector<coordinate_edge>& edges, coordinate& point, std::vector<std::vector<double>>& ridge_vectormap, double absolute_min, double absolute_range)
{
    //returns empty vector to show that it's in the closed set already
    double normalized_start_elev = (ridge_vectormap[int(point.x)][int(point.y)]-absolute_min)/(absolute_range);
    bool all_water=false;
    bool water_start = false;
    bool water_connection = false;
    if (normalized_start_elev<.4)
    {
        all_water=true;
        water_start=true;
    }
    point_and_edges to_return;
    to_return.point = point;
    to_return.coastal = false;
    std::vector<edge_information> empty_vector;
    for (int i=0; i<edges.size(); i++)
    {
        edge_information temp;
        temp.edge_index=i;
        coordinate_edge checkin = edges[i];
        double normalized_end_elev;
        if ((checkin.p1 == point) || (checkin.p2 == point))
        {
            if (checkin.p1 == point)
            {
                normalized_end_elev = (ridge_vectormap[int(checkin.p2.x)][int(checkin.p2.y)]-absolute_min)/(absolute_range);
            }
            else
            {
                normalized_end_elev = (ridge_vectormap[int(checkin.p1.x)][int(checkin.p1.y)]-absolute_min)/(absolute_range);
                
            }
            if (normalized_end_elev > .4)
            {
                all_water=false;
                double change_in_elev = (normalized_start_elev-normalized_end_elev); //want to go upwards
                temp.edge_probability=-change_in_elev;
                to_return.edge_and_probabilities.push_back(temp);
            }
            else
            {
                water_connection=true;
            }
        }
    }
    if (!all_water && water_start && water_connection)
    {
        to_return.coastal=true;
    }
    if (all_water){to_return.edge_and_probabilities=empty_vector;}
    return to_return;
}

//TODO: Add endorheic lakes to empty areas, expand size of existing lakes that get fed into
void get_to_draw (std::unordered_map<coordinate, point_and_edges, coordinate::hash>& coastal_points,
                  std::unordered_map<coordinate, point_and_edges, coordinate::hash>& land_points,
                  std::vector<coordinate_edge>& all_edges,
                  std::unordered_set<coordinate_edge, coordinate_edge::hash>& chosen_to_draw,
                  int min_size = 4)
{
    std::unordered_set<coordinate, coordinate::hash> unavailable;
    std::vector<growing_river> river_networks;
    std::random_device rd;
    std::mt19937 eng(232);
    std::uniform_real_distribution<double> threshold (0, 1);
    for (auto point_n_edge : coastal_points)
    {
        unavailable.emplace(point_n_edge.first);
        growing_river temp;
        temp.all_points.emplace(point_n_edge.first);
        for (edge_information& first_pass : point_n_edge.second.edge_and_probabilities)
        {
            coordinate_edge& of_interest = all_edges[first_pass.edge_index];
            coordinate end_point;
            if (point_n_edge.first == of_interest.p1)
            {
                end_point=of_interest.p2;
            }
            else
            {
                end_point=of_interest.p1;
            }
            if (unavailable.count(end_point)==0 && threshold(eng)<3*first_pass.edge_probability)
            {
                chosen_to_draw.emplace(of_interest);
                temp.all_points.emplace(end_point); //unnecessary? unavailable deals with this
                temp.grow_points.push_back(end_point);
                unavailable.emplace(end_point);
            }
        }
        if (temp.all_points.size()>1)
        {
            river_networks.push_back(temp);
        }
    }
    while (true)
    {
        bool no_more_grows = true;
        for (growing_river& river : river_networks)
        {
            std::vector<coordinate> current_grow_points = river.grow_points;
            river.grow_points.clear();
            if (current_grow_points.size()>0)
            {
                no_more_grows=false;
            }
            for (coordinate& grow_point : current_grow_points)
            {
                point_and_edges& grow_point_connections = land_points[grow_point]; //maybe coastal points too? see how it works out
                for (edge_information& next_passes : grow_point_connections.edge_and_probabilities)
                {
                    coordinate_edge& of_interest = all_edges[next_passes.edge_index];
                    coordinate end_point;
                    if (grow_point == of_interest.p1)
                    {
                        end_point=of_interest.p2;
                    }
                    else
                    {
                        end_point=of_interest.p1;
                    }
                    if (unavailable.count(end_point)==0 && /*threshold(eng)<(15*next_passes.edge_probability)*/ next_passes.edge_probability>0)
                    {
                        chosen_to_draw.emplace(of_interest);
                        river.all_points.emplace(end_point); //unnecessary? unavailable deals with this
                        river.grow_points.push_back(end_point);
                        unavailable.emplace(end_point);
                    }
                }
            }
        }
        if (no_more_grows)
        {
            break;
        }
    }
    //removes the short ones
    for (growing_river& river : river_networks)
    {
        if (river.all_points.size()<min_size)
        {
            //going through land points should be enough, coastal doesn't connect to coastal
            for (coordinate point : river.all_points)
            {
                for (coordinate_edge to_draw : chosen_to_draw)
                {
                    if (to_draw.p1 == point || to_draw.p2 == point) //all points are only connected on one river, thus has to be deleted
                    {
                        chosen_to_draw.erase(to_draw);
                    }
                }
            }
        }
    }
}

//always head for the sea? find closest edge, change prob to encourage
void get_to_draw (std::vector<edge_information>& edge_indices_and_probs, std::vector<coordinate_edge>& all_edges,
                  std::unordered_set<coordinate_edge, coordinate_edge::hash>& chosen_to_draw)
{

    //1 random based, 1 deterministic based til all inner points are filled.
    if (edge_indices_and_probs.size()==0)
    {
        return;
    }
    std::sort(edge_indices_and_probs.begin(), edge_indices_and_probs.end(), compare_edge_probabilities);
    double largest_probability=edge_indices_and_probs.front().edge_probability;
    int stored_sorted_index_cut_off=edge_indices_and_probs.size();
    for (int i=0; i<edge_indices_and_probs.size(); i++)
    {
        if (edge_indices_and_probs[i].edge_probability<largest_probability)
        {
            stored_sorted_index_cut_off=i;
            break;
        }
    }
    std::vector<int> edge_indices_plausible;
    for (int i=0; i<stored_sorted_index_cut_off; i++)
    {
        if (chosen_to_draw.count(all_edges[edge_indices_and_probs[i].edge_index]) == 0)
            //don't want two points to only be connected to each other
        {
            edge_indices_plausible.push_back(edge_indices_and_probs[i].edge_index);
        }
    }
    while (edge_indices_plausible.size()==0)
    {
        largest_probability=edge_indices_and_probs[stored_sorted_index_cut_off].edge_probability;
        for (int i=stored_sorted_index_cut_off; i<edge_indices_and_probs.size(); i++)
        {
            if (edge_indices_and_probs[i].edge_probability<largest_probability)
            {
                stored_sorted_index_cut_off=i;
                break;
            }
        }
        for (int i=0; i<stored_sorted_index_cut_off; i++)
        {
            if (chosen_to_draw.count(all_edges[edge_indices_and_probs[i].edge_index]) == 0)
            {
                edge_indices_plausible.push_back(edge_indices_and_probs[i].edge_index);
            }
        }
        if (largest_probability == 0 || largest_probability==edge_indices_and_probs.back().edge_probability || stored_sorted_index_cut_off == edge_indices_and_probs.size())
        {
            break;
        }
    }
    int how_many_shuffles = std::min(int(edge_indices_plausible.size()), 1);
    for (int i=0; i<how_many_shuffles; i++)
    {
        std::random_shuffle(edge_indices_plausible.begin(), edge_indices_plausible.end());
        chosen_to_draw.emplace(all_edges[edge_indices_plausible.back()]);
        edge_indices_plausible.pop_back();
    }
}


#endif /* river_creation_hpp */
