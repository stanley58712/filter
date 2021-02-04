#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define pi 3.14159265358
struct header{
	char ChunkID[4];			
	int ChunkSize;				 
	char Format[4]; 			 
	char SubChunk1ID[4];		
	int SubChunk1Size;			 
	short AudioFormat;			
	short NumChannels;			
	int SampleRate;				
	int ByteRate;				
	short BlockAlign;			
	short BitsPerSample;		
	char SubChunk2ID[4];		
	int SubChunk2Size;		
}Wave;

int main(/*int argc, char *argv[]*/){
	/*char* WavIn_Name = argv[1];//input.wav 
	char* WavOut_Name = argv[2];//output.wav
	FILE*fp_input;
	fp_input = fopen(WavIn_Name,"rb");
	FILE*fp_output;
	fp_output = fopen(WavOut_Name,"wb");*/
	
	FILE*fp_input;
	fp_input = fopen("Ascience-Fast-Piano-Add-Tones.wav","rb");
	FILE*fp_output;
	fp_output = fopen("output.wav","wb");
	
	fread(&Wave, 44, 1, fp_input);
	
	int i,j,fs,Num,Time;
	float T;
	fs = Wave.SampleRate;
	Num = Wave.SubChunk2Size/Wave.BlockAlign; //fs*T
	
	T = (float)Num/fs;//T=44.516
	Time = fs*(T+0.1);
	
	short *wave1;
	wave1 = calloc(Time,sizeof(short));
	short *wave2;
	wave2 = calloc(Time,sizeof(short));
//stereo
	for(i=0; i<Num*2; i+=2){
    	fread(&wave1[i/2],2,1,fp_input);
    	fread(&wave2[i/2],2,1,fp_input);
	}
	
	int N,P,m,M,n,k;
	P= 0.1*fs;// P = 4800
	N = 0.1*fs;// N = 4800
	M = 0.05*fs;// M = 2400
	m = (1/0.05)*T;// m = 892
	
	for(i=Num-1;i<Time;i++){
		wave1[i]=0;
		wave2[i]=0;
	}
	
	float *R1 = calloc(N,sizeof(float));
    float *I1 = calloc(N,sizeof(float));
    float *R2 = calloc(N,sizeof(float));
    float *I2 = calloc(N,sizeof(float));
	float Xn;
	int N1, N2;
	long long int LPF1,LPF2;
	float temp1,temp2;
	short *out1 = calloc(Time,sizeof(short));
	short *out2 = calloc(Time,sizeof(short));
	fwrite(&Wave,44,1,fp_output); 
//creat sine and cosine table
	float *cosine = calloc(N,sizeof(float));
	float *sine = calloc(N,sizeof(float));
	for(i=0;i<N;i++){
		cosine[i] = cos((2*pi/N)*i);
		sine[i] = sin((2*pi/N)*i);
	}
//design window types
	float *w = calloc(N,sizeof(float));;
	float *X1 = calloc(N,sizeof(float));;
	float *X2 = calloc(N,sizeof(float));;
	for(n=0;n<P;n++){
		w[n]=0.5-0.5*cos(2*pi*n/(P-1));
	}

//LPF
	LPF1 = 15000;
    N1 = LPF1*N/fs;// N1 = 1500 
    N2 = N - N1;// N2 = 4800-1500 = 3300
    printf("%d\n",m);
	int count=0;
	
	for(i=0;i<m;i++){
		count++;
		printf("%d\n",count);
	//framing and taking window
		for(k=0;k<N;k++){
			X1[k] = wave1[i*M+k]*w[k];	
			X2[k] = wave2[i*M+k]*w[k];
		}
	//DFT	
		for(k=0;k<N;k++){
			R1[k]=0;
			I1[k]=0;
			R2[k]=0;
			I2[k]=0;
			for(n=0;n<N;n++){
				j = (n * k ) % N;
				R1[k] += cosine[j]*X1[n];
				I1[k] -= sine[j]*X1[n];	
				R2[k] += cosine[j]*X2[n];
				I2[k] -= sine[j]*X2[n];
			}
		}
	//LPF
   		for(k=0; k<N; k++){
    		if(k>=N1 && k<=N2){ 
    			//R1[k] = 0.0;
    			//I1[k] = 0.0;
    			R2[k] = 0.0;
    			I2[k] = 0.0;
			}
		}
	//IDFT 
		for(k=0; k<N; k++){
			for(n=0; n<N; n++){
				j = (n * k) % N;
				temp1 += R1[n]*cosine[j] - I1[n]*sine[j];
				temp2 += R2[n]*cosine[j] - I2[n]*sine[j];
			}
			out1[i*M+k] += temp1/N;
			out2[i*M+k] += temp2/N;
			temp1 = 0.0;
			temp2 = 0.0;
		}
	}
//delete the noise below 100Hz
	int noisesize = fs*0.05;// noisesize = 2400
	short *noise = calloc(fs*0.05,sizeof(short));
	int a = fs*0.1;
	for(i=a;i<=fs*0.15;i++)
		noise[i-a] = out1[i];
	for(i=a;i<Num;i++){
		j = i % noisesize;
		out1[i] -= noise[j];
		//out2[i] -= noise[j]; 
	}
	for(i=0;i<a;i++){
		out1[i] = 0;
		//out2[i] = 0; 
	}
//output of the wave
	for(i=0;i<Num;i++){
		fwrite(&out1[i],2,1,fp_output);
		fwrite(&out2[i],2,1,fp_output);
	}
    fclose(fp_input);
    fclose(fp_output);
	printf("Program Ended\n");//µ{¦¡µ²§ô
}
	
			

