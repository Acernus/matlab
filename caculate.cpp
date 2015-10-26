#include "mex.h"
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <math.h>
#include <map>
#include <stdlib.h>
#define pi 3.1415926535
#define R0 400 //射线源到中心距离
#define R1 0 //探测器到中心距离
#define imgWidth 512 //重建图像宽度
#define imgHeight 512//重建图像高度
#define M 360 //角度
#define N 600 //探测器个数
#define iterativeTime 1 //迭代次数
#define littledelta 10
#define belta 10
using namespace std;

void initImageArray(vector<vector<double> > &img, double *data, int IM, int IN) {
	for (int i = 0; i < IM; ++i) {
		for (int j = 0; j < IN; ++j) {
			img[i][j] = data[j * IM + i];
		}
	}
}


double caculateUX(double x1, double x2) {
	if(abs(x1 - x2) < littledelta) {
		return x1 - x2;
	}
	else {
		return littledelta * x1; 
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int imgR = static_cast<int>(mxGetM(prhs[0]));
	int imgC = static_cast<int>(mxGetN(prhs[0]));
	int proR = static_cast<int>(mxGetM(prhs[1]));
	int proC = static_cast<int>(mxGetN(prhs[1]));

	if(imgR != imgHeight || imgC != imgWidth)
		return;

	double *output;
    plhs[0] = mxCreateDoubleMatrix(imgR, imgC, mxREAL);
    output = mxGetPr(plhs[0]);

	vector<vector<double> > img(imgR, vector<double>(imgC));
	vector<vector<double> > projection(proR, vector<double>(proC));

	initImageArray(img, mxGetPr(prhs[0]), imgR, imgC);// initial the array
	//获得投影数据
	initImageArray(projection, mxGetPr(prhs[1]), proR, proC);

	mexPrintf("img width: %d, height: %d \n", imgC, imgR);

	mexPrintf("projection width: %d, height: %d \n", proC, proR);

	double sintable[M], costable[M]; //每个角度的cos值

	//计算每个角度的sin和cos值
	for (int i = 0; i < M; ++i) {
		sintable[i] = sin((i + 1) * pi * 2 / M);
		costable[i] = cos((i + 1) * pi * 2 / M);
	}

	//begin the itertive loop
	//every time
	for (int i = 0; i < iterativeTime; ++i) {
		//every angle
		for (int j = 0; j < M; ++j) {
			mexPrintf("angle loop : %d\n", j + 1);
			vector<unordered_map<string, double> > lineTrack(N);//保存当前角度下射线每条射线穿过的像素点和对应的线段长
			vector<double> lineSum(N);
			vector<vector<double> > pixelSum(imgR, vector<double>(imgC));
			vector<double> delta(N);
			vector<vector<double> > deltaMutilpyPixel(imgR, vector<double>(imgC));
			for (int k = 0; k < N; ++k) {
				//计算射线源的坐标
				double x0 = R0 * costable[j], y0 = R0 * sintable[j];
				//计算探测器的射线坐标
				double x1 = -(R1 * costable[j] + (N / 2 - k) * sintable[j]), y1 = (N / 2 - k) * costable[j] - R1 * sintable[j];
				//计算斜率处理极端情况
				//视为与X轴平行
				if (abs(y1 - y0) < 1e-6) {
					if(y0 > -imgHeight / 2 && y0 < imgHeight / 2) {
						int tmpx = static_cast<int>(imgHeight / 2 - ceil(y0));
						if(abs(y0 - ceil(y0)) > 1e-6) {
							for (int tmp = 0; tmp < imgWidth; ++tmp) {
								//存入hash
								char buff[200];
								sprintf(buff, "%d %d", tmpx, tmp);
								string str(buff);
								lineTrack[j][str] = 1;
								lineSum[k] += img[tmp][tmpx];
							}
						}
					}
					//如果只有一个交点或没有，那么直接到下一条射线
					else {
						continue;
					}
				}
				else if(abs(x1 - x0) < 1e-6) {
					if(x0 > -imgWidth / 2 && x0 < imgWidth / 2) {
						int tmpy = static_cast<int>(imgWidth / 2 + floor(x0));
						if(abs(x0 - floor(x0)) > 1e-6) {
							for (int tmp = 0; tmp < imgHeight; ++tmp) {
								//存入hash
								char buff[200];
								sprintf(buff, "%d %d", tmp, tmpy);
								string str(buff);
								lineTrack[j][str] = 1;
								lineSum[k] += img[tmpy][tmp];
							}
						}
					}
					//如果只有一个交点或没有，那么直接到下一条射线
					else {
						continue;
					}
				}
				//处理一般情况
				else {
					//计算射线与重建图像的交点
					double ix0 = -imgWidth / 2, iy0 = imgHeight / 2;
                    double ix1 = imgWidth / 2, iy1 = imgHeight / 2;
                    double ix2 = -imgWidth / 2, iy2 = -imgHeight / 2;
                    double ix3 = imgWidth / 2, iy3 = -imgHeight / 2;
                    double k1 = (y1 - y0) / (x1 - x0);
                    double b = y1 - k1 * x1;
                    vector<pair<double, double> > join;
                    if (ix0 * k1 + b <= iy0 && ix0 * k1 + b >= iy2)
                        join.push_back(make_pair(ix0, k1 * ix0 + b));
                    if (ix1 * k1 + b <= iy1 && ix1 * k1 + b >= iy3)
                        join.push_back(make_pair(ix1, k1 * ix1 + b));
                    if ((iy0 - b) / k1 >= ix0 && (iy0 - b) / k1 <= ix1)
                        join.push_back(make_pair((iy0 - b) / k1, iy0));
                    if ((iy2 - b) / k1 >= ix2 && (iy2 - b) / k1 <= ix3)
                    	join.push_back(make_pair((iy2 - b) / k1, iy2));
                    //如果只有一个交点或没有，那么直接到下一条射线
                    if(join.size() < 2)
                    	continue;

                    //找到两个交点，开始遍历
                    if(join[0].first > join[1].first)
                    	std::swap(join[0], join[1]);
                    double wx0 = join[0].first, wy0 = join[0].second;
                    double wx1 = join[1].first, wy1 = join[1].second;
                    //从左到右开始遍历
                    map<double, double> m;
                    if(wx0 != ceil(wx0)) {
                    	m[wx0] = wy0;
                    }
                    for(int i = ceil(wx0); i < wx1; ++i) {
                   		if(abs(k1 * i + b - floor(k1 * i + b)) < 1e-6) {
                   			m[i] = floor(k1 * i + b);
                    	}
                    	else if(abs(k1 * i + b - ceil(k1 * i + b)) < 1e-6) {
                    		m[i] = ceil(k1 * i + b);
                    	}
                    	else
		                	m[i] = k1 * i + b;
                    }
                    if(k1 > 0) {
	                    for(int i = ceil(wy0) + 1; i < wy1; ++i) {
	                    	if(abs((i - b) / k1 - floor((i - b) / k1)) < 1e-6){
	                    		m[floor((i - b) / k1)] = i;
	                    	}
			                else if(abs((i - b) / k1 - ceil((i - b) / k1)) < 1e-6){
	                    		m[ceil((i - b) / k1)] = i;
	                    	}
	                    	else
	                    		m[(i - b) / k1] = i;
	                    }
	                }
	                else {
	                	for(int i = ceil(wy0) - 1; i > wy1; --i) {
		                    if(((i - b) / k1 - floor((i - b) / k1)) < 1e-6){
	                    		m[floor((i - b) / k1)] = i;
	                    	}
			                else if(((i - b) / k1 - ceil((i - b) / k1)) < 1e-6){
	                    		m[ceil((i - b) / k1)] = i;
	                    	}
	                    	else
	                    		m[(i - b) / k1] = i;
		                }
	                }
	                m[wx1] = wy1;

	                map<double, double>::iterator itrPrev = m.begin(), itrNext = m.begin();
	                ++itrNext;

					for(; itrNext != m.end(); ++itrPrev, ++itrNext) {
						int pixelx = floor((itrPrev->first + itrNext->first) / 2 + imgWidth / 2);
						int pixely = abs(floor((itrPrev->second + itrNext->second) / 2 - imgHeight / 2));
						//计算当前像素的位置，记录线段的长度
						char buff[200];
						sprintf(buff, "%d %d", pixelx, pixely);
						string str(buff);
								
						lineTrack[k][str] = sqrt((itrNext->first - itrPrev->first) * (itrNext->first - itrPrev->first)  + (itrNext->second - itrPrev->second) * (itrNext->second - itrPrev->second));
						lineSum[k] += lineTrack[k][str] * img[pixelx][pixely];
					}
				}

				delta[k] = projection[j][k] / lineSum[k];
			}
			mexPrintf("lineSum end\n");


			for(int k = 0; k < N; ++k) {
				for(unordered_map<string, double>::iterator itr = lineTrack[k].begin(); lineTrack[k].size() != 0 && itr != lineTrack[k].end(); ++itr) {
					string::size_type pos = itr->first.find(" ");
					string str = string(itr->first.begin(), itr->first.begin() + pos);
					int row = atoi(string(itr->first.begin(), itr->first.begin() + pos).c_str());
					int col = atoi(string(itr->first.begin() + pos + 1, itr->first.end()).c_str());
					pixelSum[row][col] += itr->second;
					deltaMutilpyPixel[row][col] += itr->second * delta[k];
				}
			}
			mexPrintf("pixelSum end\n");
			for(int k = 0; k < N; ++k) {
				for(unordered_map<string, double>::iterator itr = lineTrack[k].begin(); lineTrack[k].size() != 0 && itr != lineTrack[k].end(); ++itr) {
					string::size_type pos = itr->first.find(" ");
					int row = atoi(string(itr->first.begin(), itr->first.begin() + pos).c_str());
					int col = atoi(string(itr->first.begin() + pos + 1, itr->first.end()).c_str());

					double ux = 0;
					if(row - 1 >= 0) {
						ux += caculateUX(img[row][col], img[row - 1][col]);
					}
					if(row + 1 < imgHeight) {
						ux += caculateUX(img[row][col], img[row + 1][col]);
					}
					if(col - 1 >= 0) {
						ux += caculateUX(img[row][col], img[row][col - 1]);
					}
					if(col + 1 < imgWidth) {
						ux += caculateUX(img[row][col], img[row][col + 1]);
					}
					double c = deltaMutilpyPixel[row][col] / (pixelSum[row][col] + belta * ux);
					img[row][col] *= c;
				}
			}
			mexPrintf("loop end!\n");
		}
	}

	int k = 0;
	for(int i = 0; i < imgWidth; ++i) {
		for(int j = 0; j < imgHeight; ++j) {
			output[k++] = img[j][i];
		}
	}

}