#include "bmpBlackWhite.h"
#include "mpi.h"

/** Show log messages */
#define SHOW_LOG_MESSAGES 1

/** Enable output for filtering information */
#define DEBUG_FILTERING 0

/** Show information of input and output bitmap headers */
#define SHOW_BMP_HEADERS 0


int main(int argc, char** argv){

	tBitmapFileHeader imgFileHeaderInput;			/** BMP file header for input image */
	tBitmapInfoHeader imgInfoHeaderInput;			/** BMP info header for input image */
	tBitmapFileHeader imgFileHeaderOutput;			/** BMP file header for output image */
	tBitmapInfoHeader imgInfoHeaderOutput;			/** BMP info header for output image */
	char* sourceFileName;							/** Name of input image file */
	char* destinationFileName;						/** Name of output image file */
	int inputFile, outputFile;						/** File descriptors */
	unsigned char *outputBuffer;					/** Output buffer for filtered pixels */
	unsigned char *inputBuffer;						/** Input buffer to allocate original pixels */
	unsigned char *auxPtr;							/** Auxiliary pointer */
	unsigned int rowSize;							/** Number of pixels per row */
	unsigned int rowsBlocks;						/** Number of blocks */
	unsigned int rowsPerProcess;					/** Number of rows to be processed (at most) by each worker */
	unsigned int extraRows;
	unsigned int rowsSentToWorker;					/** Number of rows to be sent to a worker process */
	unsigned int receivedRows;						/** Total number of received rows */
	unsigned int threshold;							/** Threshold */
	unsigned int currentRow;						/** Current row being processed */
	unsigned int currentPixel;						/** Current pixel being processed */
	unsigned int outputPixel;						/** Output pixel */
	unsigned int readBytes;							/** Number of bytes read from input file */
	unsigned int writeBytes;						/** Number of bytes written to output file */
	unsigned int totalBytes;						/** Total number of bytes to send/receive a message */
	unsigned int numPixels;							/** Number of neighbour pixels (including current pixel) */
	unsigned int currentWorker;						/** Current worker process */
	unsigned int *processIDs;
	tPixelVector vector;							/** Vector of neighbour pixels */
	int imageDimensions[2];							/** Dimensions of input image */
	double timeStart, timeEnd;						/** Time stamps to calculate the filtering time */
	int size, rank, tag;							/** Number of process, rank and tag */
	MPI_Status status;								/** Status information for received messages */
	unsigned int end_process = 0;					/** Signals for end and start */
	unsigned int start_process = 1;
	unsigned char *ptrIn, *ptrOut;

		// Init
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		tag = 1;
		srand(time(NULL));

		// Check the number of processes
		if (size<2){

			if (rank == 0)
				printf ("This program must be launched with (at least) 2 processes\n");

			MPI_Finalize();
			exit(0);
		}

		// Check arguments
		if (argc != 5){

			if (rank == 0)
				printf ("Usage: ./bmpFilterDynamic sourceFile destinationFile threshold numRows\n");

			MPI_Finalize();
			exit(0);
		}

		// Get input arguments...
		sourceFileName = argv[1];
		destinationFileName = argv[2];
		threshold = atoi(argv[3]);
		rowsPerProcess = atoi(argv[4]);

		// Allocate memory for process IDs vector
		processIDs = (unsigned int *) malloc (size*sizeof(unsigned int));

		// Master process
		if (rank == 0){

			// Process starts
			timeStart = MPI_Wtime();

			// Read headers from input file
			readHeaders (sourceFileName, &imgFileHeaderInput, &imgInfoHeaderInput);
			readHeaders (sourceFileName, &imgFileHeaderOutput, &imgInfoHeaderOutput);

			// Write header to the output file
			writeHeaders (destinationFileName, &imgFileHeaderOutput, &imgInfoHeaderOutput);

			// Calculate row size for input and output images
			rowSize = (((imgInfoHeaderInput.biBitCount*imgInfoHeaderInput.biWidth) + 31) / 32 ) * 4;
			rowsBlocks = imgInfoHeaderInput.biHeight/rowsPerProcess;
			extraRows = imgInfoHeaderInput.biHeight%rowsPerProcess;

			// Show info before processing
			if (SHOW_LOG_MESSAGES){
				printf ("[MASTER] Applying filter to image %s (%dx%d) with threshold %d. Generating image %s\n", sourceFileName,
						rowSize, imgInfoHeaderInput.biHeight, threshold, destinationFileName);				
				printf ("Number of workers:%d -> Each worker calculates (at most) %d rows and %d blocks\n", size-1, rowsPerProcess, rowsBlocks);
			}

			// Show headers...
			if (SHOW_BMP_HEADERS){
				printf ("Source BMP headers:\n");
				printBitmapHeaders (&imgFileHeaderInput, &imgInfoHeaderInput);
				printf ("Destination BMP headers:\n");
				printBitmapHeaders (&imgFileHeaderOutput, &imgInfoHeaderOutput);
			}

			// Open source image
			if((inputFile = open(sourceFileName, O_RDONLY)) < 0){
				printf("ERROR: Source file cannot be opened: %s\n", sourceFileName);
				exit(1);
			}

			// Open target image
			if((outputFile = open(destinationFileName, O_WRONLY | O_APPEND, 0777)) < 0){
				printf("ERROR: Target file cannot be open to append data: %s\n", destinationFileName);
				exit(1);
			}

			// Allocate memory to copy the bytes between the header and the image data
			outputBuffer = (unsigned char*) malloc ((imgFileHeaderInput.bfOffBits-BIMAP_HEADERS_SIZE) * sizeof(unsigned char));

			// Copy bytes between headers and pixels
			lseek (inputFile, BIMAP_HEADERS_SIZE, SEEK_SET);
			read (inputFile, outputBuffer, imgFileHeaderInput.bfOffBits-BIMAP_HEADERS_SIZE);
			write (outputFile, outputBuffer, imgFileHeaderInput.bfOffBits-BIMAP_HEADERS_SIZE);

			// Allocate memory for buffers
			inputBuffer = (unsigned char*) malloc(imgInfoHeaderInput.biHeight*rowSize*sizeof(unsigned char));
			outputBuffer = (unsigned char*) realloc(outputBuffer, imgInfoHeaderInput.biHeight*rowSize*sizeof(unsigned char));
			
			// Read whole data to input buffer
			if ((readBytes = read(inputFile, inputBuffer, imgInfoHeaderInput.biHeight*rowSize)) != imgInfoHeaderInput.biHeight*rowSize) {
				showError ("Error while reading from source file");
			}
			
			ptrIn = inputBuffer;
			ptrOut = outputBuffer;
			receivedRows = 0;
			rowsSentToWorker = 0;
			currentRow = 0;

			MPI_Bcast(&imgInfoHeaderInput.biWidth, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
			MPI_Bcast(&rowSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

			// Send first rows blocks to all workers
			for (currentWorker=1; currentWorker<size && currentWorker<rowsBlocks; currentWorker++) {

				MPI_Send(&start_process, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
				MPI_Send(&currentRow, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
				MPI_Send(&rowsPerProcess, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
				MPI_Send(ptrIn, rowsPerProcess*rowSize, MPI_CHAR, currentWorker, tag, MPI_COMM_WORLD);
				printf("Send row %d and %d rows\n", currentRow, rowsPerProcess);
				ptrIn += rowsPerProcess*rowSize; // Move to next block
				rowsSentToWorker++;
				currentRow += rowsPerProcess;

			}

			// End no needed workers
			for (currentWorker=rowsBlocks; currentWorker<size; currentWorker++) {
				MPI_Send(&end_process, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
			}

			unsigned int grain, nRow;

			// Send & receive left rows
			while (receivedRows < rowsBlocks) {

				// Receive processed rows from any worker
				MPI_Recv(&nRow, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // First Row
				currentWorker = status.MPI_SOURCE;
				MPI_Recv(&grain, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD, &status); // Number of rows
				printf("Receive row %d and %d rows\n", nRow, grain);
				ptrOut = outputBuffer; // Restore pointer
				ptrOut += nRow*rowSize; // Move to current block
				MPI_Recv(ptrOut, grain*rowSize, MPI_CHAR, currentWorker, tag, MPI_COMM_WORLD, &status);
				receivedRows++;



				// Send left rows
				if (rowsSentToWorker < rowsBlocks) {
					// Send info to worker
					MPI_Send(&start_process, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(&currentRow, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(&rowsPerProcess, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(ptrIn, rowsPerProcess*rowSize, MPI_CHAR, currentWorker, tag, MPI_COMM_WORLD);
					printf("Send row %d and %d rows\n", currentRow, rowsPerProcess);
					ptrIn += rowsPerProcess*rowSize; // Move to next block
					rowsSentToWorker++;
					currentRow += rowsPerProcess;
				}

				// Check for extra rows
				else if (rowsSentToWorker == rowsBlocks && extraRows != 0) {
					// Send info to worker
					MPI_Send(&start_process, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(&currentRow, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(&extraRows, 1, MPI_UNSIGNED, currentWorker, tag, MPI_COMM_WORLD);
					MPI_Send(ptrIn, extraRows*rowSize, MPI_CHAR, currentWorker, tag, MPI_COMM_WORLD);
					printf("Send extra row %d and %d rows\n", currentRow, extraRows);
					rowsSentToWorker++;
					receivedRows--;
				}

				else MPI_Send(&end_process, 1, MPI_INT, currentWorker, tag, MPI_COMM_WORLD);

			}

			// End process
			/*for (currentWorker=1; currentWorker<size; currentWorker++) {
				MPI_Send(&end_process, 1, MPI_INT, currentWorker, tag, MPI_COMM_WORLD);
			}*/

			// Write to output file
			if ((writeBytes = write(outputFile, outputBuffer, imgInfoHeaderInput.biHeight*rowSize)) != imgInfoHeaderInput.biHeight*rowSize) {
				showError ("Error while writing to destination file");
			}

			// Close files
			close (inputFile);
			close (outputFile);

			// Process ends
			timeEnd = MPI_Wtime();

			// Show processing time
			printf("Filtering time: %f\n",timeEnd-timeStart);
		}


		// Worker process
		else {

			unsigned int width, grain;
			unsigned int signal;

			// Receive common data from master
			MPI_Bcast(&width, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
			MPI_Bcast(&rowSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
			
			// Allocate memory for buffers
			inputBuffer = (unsigned char*) malloc(rowsPerProcess*rowSize*sizeof(unsigned char));
			outputBuffer = (unsigned char*) malloc(rowsPerProcess*rowSize*sizeof(unsigned char));

			// Receive signal from master
			MPI_Recv(&signal, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			while (signal != end_process) {

				MPI_Recv(&receivedRows, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&grain, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(inputBuffer, grain*rowSize, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				// Process rows
				currentRow=0;
				rowsSentToWorker=0;

				// Process the row
				while (rowsSentToWorker < grain) {

					for (currentPixel=0; currentPixel<rowSize; currentPixel++) {

						// Current pixel
						numPixels = 0;
						vector[numPixels] = inputBuffer[currentRow+currentPixel];
						numPixels++;

						// Not the first pixel
						if (currentPixel > 0) {
							vector[numPixels] = inputBuffer[currentRow+currentPixel-1];
							numPixels++;
						}

						// Not the last pixel
						if (currentPixel < (width-1)) {
							vector[numPixels] = inputBuffer[currentRow+currentPixel+1];
							numPixels++;
						}

						outputBuffer[currentRow+currentPixel] = calculatePixelValue(vector, numPixels, threshold, DEBUG_FILTERING);						
					}

					currentRow+=rowSize;
					rowsSentToWorker++;
				}

				// Send data to master
				MPI_Send(&receivedRows, 1, MPI_UNSIGNED, 0, tag, MPI_COMM_WORLD);
				MPI_Send(&grain, 1, MPI_UNSIGNED, 0, tag, MPI_COMM_WORLD);
				MPI_Send(outputBuffer, grain*rowSize, MPI_CHAR, 0, tag, MPI_COMM_WORLD);

				// Receive signal
				MPI_Recv(&signal, 1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}

		// Finish MPI environment
		MPI_Finalize();
}
