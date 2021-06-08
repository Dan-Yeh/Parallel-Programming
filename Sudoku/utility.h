#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <omp.h>

#define BOX 3
#define N (BOX*BOX)

int parseChar(char c) {
    if (c == '.') {
        return 0;
    }
    else if (c >= '0' && c <= '9') {
        return c - '0';
    }
    else if (c >= 'A' && c <= 'P') {
        return c - 'A' + 10;
    }
    else if (c >= 'a' && c <= 'P') {
        return c - 'a' + 10;
    }
    else {
        return 0;
    }
}

bool isSafe(int grid[N][N], int row, int col, int num)
{
    for (int x = 0; x < N; x++) {
        if (grid[row][x] == num) {
            return false;
        }
    }
 
    for (int x = 0; x < N; x++) {
        if (grid[x][col] == num) {
            return false;
        }
    }
 
    int startRow = row - row % BOX;
    int startCol = col - col % BOX;
   
    for (int i = 0; i < BOX; i++) {
        for (int j = 0; j < BOX; j++) {
            if (grid[i + startRow][j + startCol] == num) {
                return false;
            }
        }
    }
 
    return true;
}

void print(int arr[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++) {
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

std::string readInput() {
    std::cout << "READY" << std::endl;

    std::string line;
    std::cin >> line;

    return line;
}
