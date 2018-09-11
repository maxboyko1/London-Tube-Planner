# London-Tube-Planner
A program for planning trips on the City of London's public transit system, complete with expected trip times and instructions for changing train lines to achieve the quickest possible trip between two stations.

Stations accepted as valid start and end points include all stations directly served by any of the rail lines operated by Transport for London, as of May 2018. This includes all 11 London Underground subway lines, as well as the Docklands Light Railway (DLR), London Overground, TfL Rail and Tramlink light rail lines. A list of valid stations, in alphabetical order, is provided in the file `stationlist.txt` for reference.

The user input to the program should be in the following format:

```
<start_station> TO <destination_station>
<start_station> TO <destination_station>
etc.
```

A sample input file to the program, `exampleinput.txt`, is provided.

For each pair of stations provided in the input file, the program will print directions for traveling from the start station to the destination in the least expected amount of time, using data on the frequency of train arrivals for each line and station, and the expected transit times between stations.

Directions given assume that all train lines are open and operating for the duration of the commuter's trip, and the expected transit times shown assume the commuter is travelling during off-peak hours.

Tested on Mac OS X 10.13.6.
