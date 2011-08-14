#ifndef MATRIXINVERSE_H
#define MATRIXINVERSE_H
#include <vector>
#include <iostream>

std::vector< std::vector< double> >
inv(const std::vector< std::vector< double> >& matrix);



std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix);


#endif // MATRIXINVERSE_H
