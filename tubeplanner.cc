/**
 * A program for planning trips on the City of London's public transit system,
 * complete with expected trip times and instructions for changing train lines
 * to achieve the quickest possible trip between two stations.
 * 
 * Stations accepted as valid start and end points include all stations
 * directly served by any of rail lines operated by Transport for London, as
 * of May 2018. This includes the 11 London Underground subway lines, as well 
 * as the Docklands Light Railway (DLR), London Overground, TfL Rail, and 
 * Tramlink light rail lines. A list of valid stations, in alphabetical order, 
 * is provided in the file stationlist.txt for reference. 
 * 
 * The user input to the program should be in the following format.
 * 
 * <start_station> TO <destination_station>
 * <start_station> TO <destination_station>
 * etc.
 * 
 * A sample input file to the program (exampleinput.txt) is provided.
 * 
 * For each pair of stations provided in the input file, the program will print
 * directions for traveling from the start station to the destination in the
 * least expected amount of time, using data on the frequency of train arrivals
 * for each line and station, and the expected transit time between adjacent 
 * stations.
 * 
 * Directions given assume that all train lines are open and operating for the 
 * duration of the commuter's trip, and the expected transit times shown assume
 * the commuter is travelling during off-peak hours.
 */

#include "tubeplanner.h"

int main ()
{
    /* List of stations we'll use to represent the graph */
    std::vector<StationVertex> station_vertices;

    /* Map of station names to another (nested) map of train lines the station 
     * is on to the corresponding StationVertex ID */
    StationLineIDsMap station_line_ids;

    read_in_train_connections(station_vertices, station_line_ids);
    read_in_interchanges(station_vertices, station_line_ids);

    int query_count = 0;
    bool errors_occurred = false;
    
    std::string query;
    std::getline(std::cin, query);

    while (!std::cin.eof() && !query.empty()) {
        std::stringstream ss(query);
        
        std::string from_station = "";
        std::string token;

        /* Read in commuter's desired starting station. */
        while (ss >> token && token.compare("TO"))
            from_station += " " + token;
        from_station.erase(0, 1);

        if (token.compare("TO")) {
            std::cerr << "ERROR: Input format must be <start_station> TO <destination_station>" 
                << std::endl;
            errors_occurred = true;

        } else {
            bool from_station_valid = station_line_ids.find(from_station) 
                != station_line_ids.end();

            if (!from_station_valid) {
                std::cerr << "ERROR: Invalid starting station: " << from_station 
                    << std::endl;
                errors_occurred = true;
            }

            std::string to_station = "";

            /* Read in commuter's desired destination station. */
            while (ss >> token)
                to_station += " " + token;
            to_station.erase(0, 1);

            bool to_station_valid = station_line_ids.find(to_station) 
                != station_line_ids.end();

            if (!to_station_valid) {
                std::cerr << "ERROR: Invalid destination station: " << to_station
                    << std::endl;
                errors_occurred = true;
            }

            if (from_station_valid && to_station_valid) {
                std::deque<int> station_path = get_fastest_route(station_vertices,
                    station_line_ids, from_station, to_station);

                std::cout << "QUERY " << ++query_count << ": " << from_station 
                    << " TO " << to_station << std::endl;
                print_directions(station_vertices, station_path);
            }
        }

        std::cout << std::endl << std::endl;
        std::getline(std::cin, query);
    }

    return errors_occurred ? 1 : 0;
}

/* 
 * Read data from global list of train connections between stations and add
 * the corresponding station edges and vertices to the graph. Store the IDs
 * for all stations (and each line that station is on) in the map for later
 * lookup.
 */
void read_in_train_connections (std::vector<StationVertex>& station_vertices, 
        StationLineIDsMap& station_line_ids)
{
    /* Initialize the list of stations with a dummy placeholder vertex to allow 
     * "real" station IDs to start at 1, and their positions in the list to 
     * match their IDs. */    
    StationVertex cur_vertex = {0, "", "", {}};
    station_vertices.push_back(cur_vertex);
    
    int station_count = 0;

    for (TrainConnection conn : train_connections) {
        int from_station_id;
        bool from_station_found = station_line_ids.find(conn.from_station) 
            != station_line_ids.end();
        bool line_found = false;

        if (from_station_found)
            line_found = station_line_ids[conn.from_station].find(conn.line)
                != station_line_ids[conn.from_station].end();

        /* If a StationVertex representing this train line at this starting
         * station does not yet exist in the list of verticies, we need to 
         * create one and add it to the list, and also add the corresponding 
         * ID to the map. Otherwise, if the StationVertex already existed, we
         * simply retrieve the corresponding ID from the map. */
        if (!from_station_found || (from_station_found && !line_found)) {
            from_station_id = ++station_count;
            cur_vertex = {from_station_id, conn.from_station, conn.line, {}};
            station_vertices.push_back(cur_vertex);

            station_line_ids[conn.from_station][conn.line] = from_station_id;
        } else {
            from_station_id = station_line_ids[conn.from_station][conn.line];
        }

        /* Do the same for the destination station. */
        int to_station_id;
        bool to_station_found = station_line_ids.find(conn.to_station) 
            != station_line_ids.end();
        line_found = false;

        if (to_station_found)
            line_found = station_line_ids[conn.to_station].find(conn.line)
                != station_line_ids[conn.to_station].end();

        if (!to_station_found || (to_station_found && !line_found)) {
            to_station_id = ++station_count;
            cur_vertex = {to_station_id, conn.to_station, conn.line, {}};
            station_vertices.push_back(cur_vertex);
            
            station_line_ids[conn.to_station][conn.line] = to_station_id;
        } else {
            to_station_id = station_line_ids[conn.to_station][conn.line];
        }

        /* Add this train connection as an edge in the start station vertex's
         * adjacency list. */
        StationEdge cur_edge = {to_station_id, conn.travel_time};
        station_vertices[from_station_id].edges.push_back(cur_edge);

        /* If this connection is symmetric, add the edge's inverse to the end
         * station's adjacency list. */
        if (conn.is_symmetric) {
            cur_edge = {from_station_id, conn.travel_time};
            station_vertices[to_station_id].edges.push_back(cur_edge);
        }
    }
}

/*
 * Read data from the global lists of interchanges between train stations and
 * add the corresponding station edges to the graph, using the map to look up
 * the IDs of each station.
 */
void read_in_interchanges (std::vector<StationVertex>& station_vertices,
        const StationLineIDsMap& station_line_ids)
{
    for (Interchange ic : interchanges) {
        int from_station_id = station_line_ids.at(ic.from_station).at(ic.from_line);
        int to_station_id = station_line_ids.at(ic.to_station).at(ic.to_line);

        StationEdge cur_edge = {to_station_id, ic.transfer_time};
        station_vertices[from_station_id].edges.push_back(cur_edge);
    }
}

/*
 * Find and return the fastest route for travelling from the specified start
 * and end stations on the provided graph, represented by a list of vertices
 * representing train stations. Use the station_line_ids map to look up the 
 * IDs on the graph representing the requested stations.
 */
std::deque<int> get_fastest_route (const std::vector<StationVertex>& station_vertices,
        const StationLineIDsMap& station_line_ids, const std::string& from_station, 
        const std::string& to_station)
{
    /* Valid destinations for the station path include all vertices 
     * representing the destination station in the graph (i.e. one for each
     * train line the station is on). */
    std::vector<int> valid_destinations;
    auto iter = station_line_ids.at(to_station).begin();
    while (iter != station_line_ids.at(to_station).end()) {
        valid_destinations.push_back(iter->second);
        ++iter;
    }
    const int ndestinations = valid_destinations.size();

    int best_time = std::numeric_limits<int>::max();
    std::deque<int> best_path;
    int closest_dest_id;

    const int nstations = station_vertices.size();
    iter = station_line_ids.at(from_station).begin();

    /* Construct a Dijkstra shortest paths spanning tree starting from each
     * vertex representing the start station in the graph. */
    while (iter != station_line_ids.at(from_station).end()) {
        int start = iter->second;

        bool in_tree[nstations + 1];
        int travel_time[nstations + 1];
        int parent[nstations + 1];

        for (int i = 0; i <= nstations; ++i) {
            in_tree[i] = false;
            travel_time[i] = std::numeric_limits<int>::max();
            parent[i] = -1;
        }

        travel_time[start] = 0;
        int v = start;

        /* The vertices in the shortest paths spanning tree are the ones for 
         * which we definitively know the shortest path and travel time from 
         * the start. Continue looping until all vertices in the graph are in 
         * the tree. */
        while (!in_tree[v]) {
            in_tree[v] = true;
            int nedges = station_vertices.at(v).edges.size();

            /* Update the distances to and parent vertices of each vertex in 
             * v's adjacency list, if the discovery and addition of v to the
             * tree has improved upon our current best known travel time to 
             * it. */
            for (int i = 0; i < nedges; ++i) {
                int endpoint = station_vertices.at(v).edges.at(i).to_id;
                int weight = station_vertices.at(v).edges.at(i).time_taken;

                if (travel_time[endpoint] > (travel_time[v] + weight)) {
                    travel_time[endpoint] = travel_time[v] + weight;
                    parent[endpoint] = v;
                }
            }

            v = 1;
            int best_dist = std::numeric_limits<int>::max();

            /* Choose the next value of v to be the vertex that has the 
             * shortest travel time from the start, and is not already in the
             * tree (if one exists). */
            for (int i = 1; i < nstations; ++i) {
                if (!in_tree[i] && (best_dist > travel_time[i])) {
                    best_dist = travel_time[i];
                    v = i;
                }
            }
        }

        int cur_best_time = std::numeric_limits<int>::max();
        int cur_closest_dest_id;

        /* Find the closest valid destination from the current start vertex. */
        for (int i = 0; i < ndestinations; ++i) {
            if (cur_best_time > travel_time[valid_destinations[i]]) {
                cur_best_time = travel_time[valid_destinations[i]];
                cur_closest_dest_id = valid_destinations[i];
            }
        }

        /* If this beats our overall fastest route between the two stations,
         * update the overall fastest route with the new best time and station 
         * path. */
        if (best_time > cur_best_time) {
            best_time = cur_best_time;
            closest_dest_id = cur_closest_dest_id;

            v = closest_dest_id;
            best_path.clear();

            /* Construct the new best path by following the tree's parent 
             * pointers from the destination vertex back to the start. */
            while (v != start) {
                best_path.push_front(v);
                v = parent[v];
            }

            best_path.push_front(v);
        }

        ++iter;
    }

    return best_path;
}

/*
 * Given a list of station IDs representing a path through the graph, print
 * user-readable directions for travelling from the beginning to the end of
 * the path, using the graph (i.e. list of station vertices) to look up the
 * name and train line of each station on the path as necessary, as well as
 * to calculate travel times between stations.
 */
void print_directions (const std::vector<StationVertex>& station_vertices,
        const std::deque<int>& station_path)
{
    const int path_len = station_path.size();
    
    int pos = 0;
    int total_travel_time = 0;
    int step_count = 0;
    int travel_time;

    /* Loop until we have arrived at the final station. */
    while (pos < path_len - 1) {

        /* If this is not the first station on the path, then it must be the 
         * endpoint of an interchange starting at the station immediately 
         * before it on the station path, and we must add the time this 
         * interchange takes to the total travel time. */
        if (pos) {
            travel_time = get_travel_time(station_vertices, station_path.at(pos - 1),
                station_path.at(pos));

            if (travel_time < 0) {
                std::cerr << "ERROR: Invalid station path passed to directions!" 
                    << std::endl;
                return;
            } else {
                total_travel_time += travel_time;
            }
        }

        std::string cur_line = station_vertices.at(station_path.at(pos)).line;
        
        bool is_initial_interchange = !pos && cur_line.compare(
                station_vertices.at(station_path.at(pos + 1)).line);
        
        bool is_final_interchange = (pos == path_len - 2) && cur_line.compare(
                station_vertices.at(station_path.at(pos + 1)).line);

        if (is_initial_interchange && !is_final_interchange) {
            std::cout << std::endl << ++step_count << ". From " 
                << station_vertices.at(station_path.at(pos)).name 
                << " station, proceed on foot to the "
                << station_vertices.at(station_path.at(pos + 1)).line
                << " line departure platform at the nearby "
                << station_vertices.at(station_path.at(pos + 1)).name 
                << " station. (" << total_travel_time << " minutes)" 
                << std::endl;
            
            ++pos;
        
        } else if (is_final_interchange) {
            std::cout << std::endl << ++step_count << ". From "
                << station_vertices.at(station_path.at(pos)).name 
                << " station, proceed on foot to the nearby "
                << station_vertices.at(station_path.at(pos + 1)).name
                << " station. (" << total_travel_time << " minutes)"
                << std::endl;
            
            travel_time = get_travel_time(station_vertices, station_path.at(pos),
                station_path.at(pos + 1));
                
            if (travel_time < 0) {
                std::cerr << "ERROR: Invalid station path passed to directions!" 
                    << std::endl;
                return;
            } else {
                total_travel_time += travel_time;
            }

            ++pos;
        } 

        if (!is_initial_interchange && !is_final_interchange) {
            std::cout << std::endl << ++step_count << ". Board the train at " 
                << station_vertices.at(station_path.at(pos)).name << " station on the " 
                << station_vertices.at(station_path.at(pos)).line << " line, ";

            int pass_stations_start = ++pos;

            /* Find the path position of the final station on the current line. */
            while ((pos < path_len - 1) && !cur_line.compare(station_vertices.at(
                station_path.at(pos + 1)).line))
                    ++pos;

            std::cout << "heading toward " << station_vertices.at(station_path.at(pos)).name
                << " station. (" << total_travel_time << " minutes)" << std::endl;

            if (pass_stations_start < pos)
                std::cout << "Passing stations:" << std::endl;
            
            /* Print the name of each passing station, and add the transit time 
             * between it and its predecessor to the total travel time. */
            for (int i = pass_stations_start; i < pos; ++i) {
                travel_time = get_travel_time(station_vertices, station_path.at(i - 1),
                    station_path.at(i));
                
                if (travel_time < 0) {
                    std::cerr << "ERROR: Invalid station path passed to directions!" 
                        << std::endl;
                    return;
                } else {
                    total_travel_time += travel_time;
                    std::cout << "-> " << station_vertices.at(station_path.at(i)).name 
                        << " (" << total_travel_time << " minutes)" << std::endl;
                }
            }

            /* Now add the travel time to the final station on the line to the total. */
            travel_time = get_travel_time(station_vertices, station_path.at(pos - 1),
                station_path.at(pos));
                
            if (travel_time < 0) {
                std::cerr << "ERROR: Invalid station path passed to directions!" 
                    << std::endl;
                return;
            } else {
                total_travel_time += travel_time;
            }

            /* If we have not yet arrived at the final station on the path, the 
             * commuter will now need to interchange on foot to a connecting 
             * station, which may itself be the final destination, or simply be
             * a place to switch to another train line. */
            if (pos < path_len - 1) {
                std::cout << std::endl << ++step_count << ". Disembark at " 
                    << station_vertices.at(station_path.at(pos)).name
                    << " station and proceed on foot to the ";

                is_final_interchange = (pos == path_len - 2) && cur_line.compare(
                        station_vertices.at(station_path.at(pos + 1)).line);

                if (is_final_interchange) {
                    std::cout << "nearby " << station_vertices.at(station_path.at(pos + 1)).name
                        << " station. (" << total_travel_time << " minutes)" << std::endl;

                    travel_time = get_travel_time(station_vertices, station_path.at(pos),
                        station_path.at(pos + 1));
                        
                    if (travel_time < 0) {
                        std::cerr << "ERROR: Invalid station path passed to directions!" 
                            << std::endl;
                        return;
                    } else {
                        total_travel_time += travel_time;
                    }

                } else {
                    std::string cur_station = station_vertices.at(station_path.at(pos)).name;
                    std::string next_station = station_vertices.at(station_path.at(pos + 1)).name;

                    /* Notify the commuter whether or not they are interchanging to a
                     * different station, or simply a different train line departing
                     * from the same station. */
                    if (cur_station.compare(next_station)) {
                        std::cout << station_vertices.at(station_path.at(pos + 1)).line
                            << " line departure platform at the nearby "
                            << station_vertices.at(station_path.at(pos + 1)).name
                            << " station. ";
                    } else {
                        std::cout << "same station's "
                            << station_vertices.at(station_path.at(pos + 1)).line
                            << " line departure platform. ";
                    }

                    std::cout << "(" << total_travel_time << " minutes)" << std::endl;                
                }

                ++pos;
            }
        }
    }

    std::cout << std::endl << ++step_count << ". Arrive at destination station " 
        << station_vertices.at(station_path.at(pos)).name << ". (" 
        << total_travel_time << " minutes)" << std::endl;
}

/* 
 * Helper function for print_directions() that returns the time it takes to
 * travel between the (directly connected) stations represented on the graph
 * by start_id and end_id. If the stations are not directly connected, return
 * an error code of -1.
 */
int get_travel_time (const std::vector<StationVertex>& station_vertices,
        int start_id, int end_id)
{
    auto iter = station_vertices.at(start_id).edges.begin();
    while (iter != station_vertices.at(start_id).edges.end()) {
        if (iter->to_id == end_id)
            return iter->time_taken;
        ++iter;
    }
    return -1;
}
