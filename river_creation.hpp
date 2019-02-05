//
//  river_creation.hpp
//  fun with noise
//
//  Created by asdfuiop on 2/5/19.
//  Copyright © 2019 asdfuiop. All rights reserved.
//

#ifndef river_creation_hpp
#define river_creation_hpp

#include "3rd party libs/delaunay-triangulation/delaunay.h"

struct delaunay_point_hasher
{
    std::size_t operator()(Vector2<double> const& e) const
    {
        return e.x*53 + e.y*31;
    }
};
struct edge_hasher
{
    std::size_t operator()(Edge<double> const& e) const
    {
        return (e.p1.x*53 + e.p1.y*31 + e.p2.x*53 + e.p2.y*31);
    }
};

//should redo this to ensure no doubles, use unordered set.
struct edge_information
{
    int edge_index;
    double edge_probability;
    bool on_coast;
};

bool compare_edge_probabilities(edge_information& i, edge_information& j)
{
    return i.edge_probability > j.edge_probability;
}

struct point_and_edges
{
    Vector2<double> point;
    std::vector<edge_information> edge_and_probabilities;
};

std::vector<int> find_edges_involving_point(std::vector<Edge<double>>& edges, Vector2<double>& point)
{
    std::vector<int> to_return_indices;
    for (int i=0; i<edges.size(); i++)
    {
        Edge<double> checkin = edges[i];
        if ((checkin.p1 == point) || (checkin.p2 == point))
        {
            to_return_indices.push_back(i);
        }
    }
    return to_return_indices;
}

std::vector<edge_information> find_edges_involving_point(std::vector<Edge<double>>& edges, Vector2<double>& point, std::vector<std::vector<double>>& ridge_vectormap, double absolute_min, double absolute_max)
{
    //returns empty vector to show that it's in the closed set already
    double normalized_start_elev = (ridge_vectormap[int(point.x)][int(point.y)]-absolute_min)/(absolute_max-absolute_min);
    bool all_water=false;
    bool water_start = false;
    if (normalized_start_elev<.4)
    {
        all_water=true;
        water_start=true;
    }
    std::vector<edge_information> to_return;
    std::vector<edge_information> empty_vector;
    for (int i=0; i<edges.size(); i++)
    {
        edge_information temp;
        temp.edge_index=i;
        Edge<double> checkin = edges[i];
        double normalized_end_elev;
        if ((checkin.p1 == point) || (checkin.p2 == point))
        {
            if (checkin.p1 == point)
            {
                normalized_end_elev = (ridge_vectormap[int(checkin.p2.x)][int(checkin.p2.y)]-absolute_min)/(absolute_max-absolute_min);
            }
            else
            {
                normalized_end_elev = (ridge_vectormap[int(checkin.p1.x)][int(checkin.p1.y)]-absolute_min)/(absolute_max-absolute_min);
                
            }
            if (normalized_end_elev > .4)
            {
                all_water=false;
                temp.on_coast=false;
                if (water_start)
                {
                    temp.on_coast = true;
                }
                double change_in_elev = normalized_start_elev-normalized_end_elev; //want to go upwards
                //std::cout << change_in_elev <<"\n";
                temp.edge_probability=change_in_elev;
                to_return.push_back(temp);
            }
        }
    }
    if (all_water){return empty_vector;}
    return to_return;
}

//always head for the sea? find closest edge, change prob to encourage
void get_to_draw (std::vector<edge_information>& edge_indices_and_probs, std::vector<Edge<double>>& all_edges,
                  std::unordered_set<Edge<double>, edge_hasher>& chosen_to_draw)
{
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
