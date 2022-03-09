#include <stdio.h>
#include <cuda.h>

typedef unsigned char u8;

typedef struct cell {
    u8 state;
    size_t* neighbor;
    u8 neighborSize;
} cell;

typedef enum type {life, koNeiman, koNeimanMur, koMur} type;

u8*** hStates;
size_t hX, hY, hZ;
type hT;

__device__ u8* dStates;
__device__ size_t *pdX, *pdY, *pdZ;
__device__ type *pdT;

__device__ cell* dCurrent;
__device__ cell* dNext;
__device__ size_t* pdFullSize;

void readInput(const char* inputFile)
{
    FILE* input = fopen(inputFile, "r");
    if (input == NULL) {
        printf("Can't open file %s\n", inputFile);
        exit(-1);
    }
    
    u8 firstLine = 1;
    const int LINE_SIZE = 100;
    char line[LINE_SIZE];
    size_t x, y, z;
    u8 cellState;
    int typeNumber;
    
    size_t i;
    size_t wordStart;
    u8 inWord;
    u8 separator;
    char* word;
    size_t wordSize;
    size_t wordCount;
    
    while (fgets(line, LINE_SIZE, input)) {
        wordCount = 0;
        wordStart = 0;
        inWord = 0;
        i = 0;
        if (firstLine) {
            while (line[i] != '\0') {
                separator = (line[i] == ' ' || line[i] == '\n') ? 1 : 0;
                if (inWord) {
                    if (separator) {
                        inWord = 0;
                        wordSize = i - wordStart;
                        word = (char*) malloc(wordSize + 1);
                        memcpy(word, line + wordStart, (i - wordStart) * sizeof(char));
                        word[wordSize] = '\0';
                        
                        switch (wordCount) {
                            case 0:
                                hX = (size_t)atoi(word);
                                break;
                            case 1:
                                hY = (size_t)atoi(word);
                                break;
                            case 2:
                                hZ = (size_t)atoi(word);
                                break;
                            case 3:
                                typeNumber = atoi(word);
                                switch (typeNumber) {
                                    case 1:
                                        hT = life;
                                        break;
                                    case 2:
                                        hT = koNeiman;
                                        break;
                                    case 3:
                                        hT = koNeimanMur;
                                        break;
                                    case 4:
                                        hT = koMur;
                                        break;
                                    default:
                                        printf("Wrong type of simulation: %d\n", typeNumber);
                                        exit(-1);
                                }
                                break;
                            default:
                                printf("Too much words in a line: %s\n", line);
                                exit(-1);
                        }
                        
                        free(word);
                        wordCount++;
                    }
                }
                else if (!separator) {
                    inWord = 1;
                    wordStart = i;
                }
                i++;
            }
            firstLine = 0;
            
            hStates = (u8***) malloc(hX * sizeof(u8**));
            for (x = 0; x < hX; x++) {
                hStates[x] = (u8**) malloc(hY * sizeof(u8*));
                for (y = 0; y < hY; y++) {
                    hStates[x][y] = (u8*) malloc(hZ * sizeof(u8));
                    for (z = 0; z < hZ; z++)
                        hStates[x][y][z] = 0;
                }
            }
        }
        else {
            while (line[i] != '\0') {
                separator = (line[i] == ' ' || line[i] == '\n') ? 1 : 0;
                if (inWord) {
                    if (separator) {
                        inWord = 0;
                        wordSize = i - wordStart;
                        word = (char*) malloc(wordSize + 1);
                        memcpy(word, line + wordStart, (i - wordStart) * sizeof(char));
                        word[wordSize] = '\0';
                        
                        switch (wordCount) {
                            case 0:
                                cellState = (u8)atoi(word);
                                if (cellState == 0)
                                    goto stop;
                                break;
                            case 1:
                                x = (size_t)atoi(word);
                                break;
                            case 2:
                                y = (size_t)atoi(word);
                                break;
                            case 3:
                                z = (size_t)atoi(word);
                                break;
                            default:
                                printf("Too much words in a line: %s\n", line);
                                exit(-1);
                        }
                        free(word);
                        wordCount++;
                    }
                }
                else if (!separator) {
                    inWord = 1;
                    wordStart = i;
                }
                i++;
            }
            hStates[x][y][z] = cellState;
        }
    }
    
stop:
    fclose(input);
}

void passStatesToDevice(u8 ***hStates, type *hT, size_t *hX, size_t *hY, size_t *hZ)
{
    size_t i, j;
    cudaMalloc((void**)&dStates, hX * hY * hZ * sizeof(u8));
    for (i = 0; i < hX; i++)
        for (j = 0; j < hY; j++)
            cudaMemcpy(&dStates[i * hY * hZ + j * hZ], hStates[i][j], hZ * sizeof(u8), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&pdX, sizeof(size_t));
    cudaMemcpy(pdX, &hX, sizeof(size_t), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&pdY, sizeof(size_t));
    cudaMemcpy(pdY, &hY, sizeof(size_t), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&pdZ, sizeof(size_t));
    cudaMemcpy(pdZ, &hZ, sizeof(size_t), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&pdT, sizeof(type));
    cudaMemcpy(pdT, &hT, sizeof(type), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&pdFullSize, sizeof(size_t));
    size_t size = hX * hY * hZ;
    cudaMemcpy(pdFullSize, &size, sizeof(size_t), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**)&dCurrent, size * sizeof(cell));
    cudaMalloc((void**)&dNext, size * sizeof(cell));
    
    u8 neighborSize = 6;
    for (i = 0; i < size; i++) {
        cudaMemcpy(&dCurrent[i].neighborSize, &neighborSize, sizeof(u8), cudaMemcpyHostToDevice);
        cudaMemcpy(&dNext[i].neighborSize, &neighborSize, sizeof(u8), cudaMemcpyHostToDevice);
                
        size_t* tmp1;
        size_t* tmp2;
        cudaMalloc((void**)&tmp1, neighborSize * sizeof(size_t));
        cudaMalloc((void**)&tmp2, neighborSize * sizeof(size_t));
        cudaMemcpy(&dCurrent[i].neighbor, &tmp1, neighborSize * sizeof(size_t), cudaMemcpyDeviceToDevice);
        cudaMemcpy(&dNext[i].neighbor, &tmp2, neighborSize * sizeof(size_t), cudaMemcpyDeviceToDevice);
    }
}

__device__ void idx3to1(size_t x, size_t y, size_t z, size_t* i)
{
    *i = x * (*pdY) * (*pdZ) + y * (*pdZ) + z;
}

__device__ void idx1to3(size_t i, size_t* x, size_t* y, size_t* z)
{
    *x = i / (*pdY) / (*pdZ);
    i -= (*x) * (*pdY) * (*pdZ);
    *y = i / (*pdZ);
    i -= (*y) * (*pdZ);
    *z = i;
}

__device__ size_t plus(size_t i, size_t* max)
{
    return (i == *max - 1) ? 0 : ++i;
}

__device__ size_t minus(size_t i, size_t* max)
{
    return (i == 0) ? *max - 1 : --i;
}

__global__ void transformStatesIntoCells()
{
    int idx = threadIdx.x;
    while (idx <= *pdFullSize) {
        dCurrent[idx].state = dNext[idx].state = dStates[idx];
        
        size_t x, y, z;
        idx1to3(idx, &x, &y, &z);
        
        size_t xn[2], yn[2], zn[2];
        xn[0] = minus(x, pdX);
        xn[1] = plus(x, pdX);
        yn[0] = minus(y, pdY);
        yn[1] = plus(y, pdY);
        zn[0] = minus(z, pdZ);
        zn[1] = plus(z, pdZ);
        
        int i, j, k;
        size_t neighborIdx;
        int neighborCount = 0;
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                for (k = 0; k < 2; k++) {
                    idx3to1(xn[i], yn[j], zn[k], &neighborIdx);
                    dCurrent[idx].neighbor[neighborCount] = dNext[idx].neighbor[neighborCount] = neighborIdx;
                    neighborCount++;
                }
        
        idx += blockDim.x;
    }
    __syncthreads();
}

__global__ void calc()
{
    int idx = threadIdx.x;
    while (idx <= *pdFullSize) {
        u8 s = 0;
        int i;
        for (i = 0; i < dCurrent[idx].neighborSize; i++)
            if (dCurrent[dCurrent[idx].neighbor[i]].state)
                s++;
        
        if (dCurrent[idx].state) {
            if (s < 4)
                dNext[idx].state = 0;
        }
        else {
            if (s >= 6)
                dNext[idx].state = 1;
        }
        idx += blockDim.x;
    }
    __syncthreads();
    
    idx = threadIdx.x;
    while (idx <= *pdFullSize) {
        dCurrent[idx].state = dNext[idx].state;
        idx += blockDim.x;
    }
    __syncthreads();
}

__global__ void transformCellsIntoStates()
{
    int idx = threadIdx.x;
    while (idx <= *pdFullSize) {
        dStates[idx] = dCurrent[idx].state;
        idx += blockDim.x;
    }
    __syncthreads();
}

void getDataFromDevice(size_t nThreads)
{
    transformCellsIntoStates<<<1, nThreads>>>();
    
    size_t i, j;
    for (i = 0; i < hX; i++)
        for (j = 0; j < hY; j++)
            cudaMemcpy(hStates[i][j], &dStates[i * hY * hZ + j * hZ], hZ * sizeof(u8), cudaMemcpyDeviceToHost);
}

void print(const char* outputFile)
{
    FILE* output = fopen(outputFile, "a");
    if (output == NULL) {
        printf("Can't open file %s\n", outputFile);
        exit(-1);
    }
    
    size_t i, j, k;
    for (i = 0; i < hX; i++)
        for (j = 0; j < hY; j++)
            for (k = 0; k < hZ; k++)
                fprintf(output, "%d %ld %ld %ld\n", hStates[i][j][k], i, j, k);
    fprintf(output, "0 0 0 0\n");
    fclose(output);
}

void printResults(const char* outputFile, size_t nThreads)
{
    getDataFromDevice(nThreads);
    print(outputFile);
}

void clean()
{
    cudaFree(dStates);
    cudaFree(pdX);
    cudaFree(pdY);
    cudaFree(pdZ);
    cudaFree(pdT);
    
    /*size_t i;
    for (i = 0; i < hX * hY * hZ; i++) {
        cudaFree(dCurrent[i].neighbor);
        cudaFree(dNext[i].neighbor);
    }*/
    cudaFree(dCurrent);
    cudaFree(dNext);
    cudaFree(pdFullSize);
    
    size_t x, y;
    for (x = 0; x < hX; x++) {
        for (y = 0; y < hY; y++)
            free(hStates[x][y]);
        free(hStates[x]);
    }
    cudaFree(hStates);
}

void gameOfLife(const char* inputFile, int nSteps, int outputInterval, const char* outputFile)
{
    readInput(inputFile);
    printf("Input file has been read\n");
    passStatesToDevice();
    printf("States have been copied to device\n");
    
    int device;
    cudaGetDevice(&device);
    
    struct cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);
    
    transformStatesIntoCells<<<1, prop.maxThreadsDim[0]>>>();
    printf("Neighbors have been set\n");
    
    FILE* output = fopen(outputFile, "w");
    if (output == NULL) {
        printf("Can't open file %s\n", outputFile);
        exit(-1);
    }
    fclose(output);
    
    print(outputFile);
    printf("Output\n");
    
    int i;
    for (i = 1; i <= nSteps; i++) {
        calc<<<1, prop.maxThreadsDim[0]>>>();
        printf("Step %d\n", i);
        if (i % outputInterval == 0) {
            printResults(outputFile, prop.maxThreadsDim[0]);
            printf("Output\n");
        }
    }
    
    clean();
    printf("Memory has been set free\n");
}

int main(int argc, const char * argv[])
{
    if (argc != 5)
        printf("Usage: %s inputFile nSteps outputInterval outputFile\n", argv[0]);
    else
    {
        const char* inputFile = argv[1];
        int nSteps = atoi(argv[2]);
        int outputInterval = atoi(argv[3]);
        const char* outputFile = argv[4];
        
        gameOfLife(inputFile, nSteps, outputInterval, outputFile);
    }
    return 0;
}
