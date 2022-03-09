#ifndef COMMON_H
#define COMMON_H

#include <vector>

typedef std::pair< std::pair<int,int> , std::pair<int,int> > coord;
typedef std::vector<coord> frame;
typedef std::vector<frame> frames;
#define W first.first
#define X first.second
#define Y second.first
#define Z second.second

#endif // COMMON_H
