#include "utility.h"
#include <omp.h>

int task_counter = 0;
bool solved = false;

bool isSafe(uint8_t grid[N][N], int row, int col, int num)
{
    for (int x = 0; x < N; x++)
    {
        if (grid[row][x] == num || grid[x][col] == num)
        {
            return false;
        }
    }

    int startRow = row - row % BOX;
    int startCol = col - col % BOX;

    for (int i = 0; i < BOX; i++)
    {
        for (int j = 0; j < BOX; j++)
        {
            if (grid[i + startRow][j + startCol] == num)
            {
                return false;
            }
        }
    }

    return true;
}

void print(uint8_t arr[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::cout << (int)arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

bool solveSuduko(uint8_t grid[N][N], int row, int col)
{
    if (solved)
    {
        return false;
    }

    if (row == N - 1 && col == N)
    {
        print(grid);
        return true;
    }

    if (col == N)
    {
        row++;
        col = 0;
    }

    if (grid[row][col] > 0)
    {
        return solveSuduko(grid, row, col + 1);
    }

    bool return_val = false;
    for (int num = 1; num <= N; num++)
    {
        if (return_val) continue;
        if (isSafe(grid, row, col, num))
        {
            grid[row][col] = num;
            if (row < 2)
            {
                uint8_t loc_grid[N][N];
                memcpy(loc_grid, grid, N * N * sizeof(uint8_t));
                #pragma omp taskgroup
                #pragma omp task shared(solved)
                {
                    if (solveSuduko(loc_grid, row, col + 1))
                    {
                        solved = true;
                        return_val = true;
                        #pragma omp cancel taskgroup
                    }
                }
            }
            else
            {
                if (solveSuduko(grid, row, col + 1))
                {
                    return true;
                }
            }
            grid[row][col] = 0;
        }
    }
    #pragma omp taskwait
    return false;
}

int main()
{
    std::cerr << omp_get_proc_bind() << std::endl;

    uint8_t board[N][N];

    // This starts the timer.
    std::string line = readInput();
    ;

    int counter = 0;
    char ch;
    std::istringstream iss(line);
    while (iss.get(ch))
    {
        board[counter / N][counter % N] = parseChar(ch);
        counter += 1;
    }

    #pragma omp parallel num_threads(32)
    {
        #pragma omp single
        solveSuduko(board, 0, 0);
    }
    // This stops the timer.
    std::cout << std::endl
              << "DONE" << std::endl;

    return 0;
}
