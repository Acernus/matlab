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
#include <fstream>
#define pi 3.1415926535
#define R0 388 //射线源到中心距离
#define R1 119 //探测器到中心距离
#define imgWidth 256 //重建图像宽度
#define imgHeight 256//重建图像高度
#define M 400 //角度
#define N 256 //探测器个数
#define iterativeTime 1 //迭代次数
#define offset 0 //水平方向校正
using namespace std;


const string filename = "E:\\ml_em_imgs\\ml_em_img_";

struct BIN_HEADER {	//********************* *.BIN file header struct
	char	s[492];		// Reserved
	double	min;		// Minimal value of data
	double	max;		// Maximal value of data
	int		width;		// Width of data
	int     height;		// Height of data
	int     depth;		// Depth of data (slices)
};
BIN_HEADER dataheader;
void save(int time, vector<vector<double> > &a) {
	FILE *fp;
	string file;
	file = filename + to_string(time) + ".bin";
	fopen_s(&fp, file.c_str(), "wb");
	if (fp == NULL) {
		return;
	}
	dataheader.min = a[0][0];
	dataheader.max = a[0][0];
	dataheader.height = a.size();
	dataheader.width = a[0].size();
	dataheader.depth = 1;
	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < a[0].size(); ++j) {
			if (dataheader.min > a[i][j]) dataheader.min = a[i][j];
			if (dataheader.max < a[i][j]) dataheader.max = a[i][j];
			fwrite(&a[i][j], 8, 1, fp);
		}	
	}
	fwrite(&dataheader, sizeof(BIN_HEADER), 1, fp);
	fclose(fp);
}


void initImageArray(vector<vector<double> > &img, double *data, int IM, int IN) {
	for (int i = 0; i < IM; ++i) {
		for (int j = 0; j < IN; ++j) {
			img[i][j] = data[j * IM + i];
		}
	}
}

bool cmp(pair<double, double> &a, pair<double, double> &b) {
	return a.first < b.first;
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
	//从1度到360度
	for (int i = 0; i < M; ++i) {
		sintable[i] = sin((i + 1) * pi * 0.9 * 2 / 360);
		costable[i] = cos((i + 1) * pi * 0.9 * 2 / 360);
	}


	//temp buff
	char buff[20];
	//hash string

	string str;

	for (int i = 0; i < iterativeTime; ++i) {
		//every angle

		vector<vector<double> > pixelSum(imgR, vector<double>(imgC));
		vector<vector<double> > deltaMutilpyPixel(imgR, vector<double>(imgC));
		for (int j = 0; j < M; ++j) {
			mexPrintf("angle loop : %d\n", j + 1);
			vector<unordered_map<string, double> > lineTrack(N);//保存当前角度下射线每条射线穿过的像素点和对应的线段长
			vector<double> lineSum(N);
			vector<double> delta(N);
			//坐标系重建图像左下角
			//计算射线源的坐标
			double x0 = R0 * costable[j] + imgWidth / 2, y0 = R0 * sintable[j] + imgHeight / 2;

			for (int k = 0; k < N; ++k) {
				//计算探测器的射线坐标
				double x1 = -(R1 * costable[j] + (N / 2 - k + offset) * sintable[j]) + imgWidth / 2, y1 = -(R1 * sintable[j] - (N / 2 - k + offset) * costable[j]) + imgHeight / 2;
				//计算斜率处理极端情况
				//视为与X轴平行
				double k1, b, xmin, xmax, ymin, ymax;

				if (abs(y1 - y0) < 1e-6) {
					//y0在图像的范围内
					if(y0 > 0 && y0 < imgHeight) {
						int tmpx = static_cast<int>(floor(y0));
						//y0不在图像的边界上
						if(abs(y0 - floor(y0)) > 1e-6) {
							for (int i = 0; i < imgWidth; ++i) {
								//存入hash
								memset(buff, '\0', sizeof(buff));
								sprintf(buff, "%d %d", tmpx, i);
								str = buff;
								lineTrack[k][str] = 1;
								lineSum[k] += img[tmpx][i] * 1;
							}
						} else {
							continue;
						}
					}
					//如果在边界，那么直接到下一条射线
					else {
						continue;
					}
				}
				else if(abs(x1 - x0) < 1e-6) {
					//x0在图像的范围内
					if(x0 > 0 && x0 < imgWidth) {
						int tmpy = static_cast<int>(floor(x0));
						if(abs(x0 - floor(x0)) > 1e-6) {
							for (int i = 0; i < imgHeight; ++i) {
								//存入hash
								memset(buff, '\0', sizeof(buff));
								sprintf(buff, "%d %d", i, tmpy);
								str = buff;
								lineTrack[k][str] = 1;
								lineSum[k] += img[i][tmpy] * 1;
							}
						} else {
							continue;
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
					k1 = (y1 - y0) / (x1 - x0);
        			b = y0 - k1 * x0;
                    if( k1 > 0) {
			            xmin = (-b / k1) < 0 ? 0 : (-b / k1) ;
			            xmax = ((imgHeight - b) / k1) > imgWidth ? imgWidth : (imgHeight - b) / k1;
			            ymin = xmin * k1 + b;
			            ymax = xmax * k1 + b;
			        }
			        if( k1 < 0) {
			            xmin = (imgHeight - b) / k1 < 0 ? 0 : (imgHeight - b) / k1 ;
			            xmax = -b / k1 > imgWidth ? imgWidth : -b / k1;
			            ymin = xmax * k1 + b;
			            ymax = xmin * k1 + b;
			        }
			        //如果射线与重建图像没有交点
			        if(xmin >= imgWidth || xmax <= 0) {
			            continue;
			        }
                    else {
			        	vector<pair<double, double> > v;
			        	//插入每个相交的x点
			            for(int i = static_cast<int>(ceil(xmin)); i <= static_cast<int>(floor(xmax)); ++i) {
							v.push_back(make_pair(i, k1 * i + b));
			            }

			            if(k1 > 0) {
			                for(int i = static_cast<int>(ceil(ymin)); i <= static_cast<int>(floor(ymax)); ++i) {
								v.push_back(make_pair((i - b) / k1, i));
			                }
			            }
			            else if(k1 < 0) {
			                for(int i = static_cast<int>(floor(ymax)); i >= static_cast<int>(ceil(ymin));--i) {
								v.push_back(make_pair((i - b) / k1, i));
			                }
			            }
			            
			            //添加第一个元素
			            //如果xmin不为ceil接近
			            if(abs(xmin - ceil(xmin)) > 1e-6) {
			            	v.push_back(make_pair(xmin, k1 * xmin + b));
			            }
			            //添加最后一个元素
			            //如果xmax不为floor接近
			            if(abs(xmax - floor(xmax)) > 1e-6) {
			            	v.push_back(make_pair(xmax, k1 * xmax + b));
			            }

			            //排序清除相同的元素。
						sort(v.begin(), v.end(), cmp);
			            vector<pair<double, double> > tmpVector;
			          	vector<pair<double, double> >::iterator tmp = v.begin();
						tmpVector.push_back(*tmp);
			            for(vector<pair<double, double> >::iterator itr = v.begin() + 1; itr != v.end(); ++itr) {
			            	if(!(abs(itr->first - tmp->first) < 1e-6 && abs(itr->second - tmp->second) < 1e-6)) {
			            		tmp = itr;
								tmpVector.push_back(*tmp);
			            	}
			            }
			            v = tmpVector;
			            
			            for(int i = 0; i < v.size() - 1; ++i) {
			            	int col = floor((v[i].first + v[i + 1].first) / 2);
							int row = imgHeight - ceil((v[i].second + v[i + 1].second) / 2);
							memset(buff, '\0', sizeof(buff));
							sprintf(buff, "%d %d", row, col);
							str = buff;

							lineTrack[k][str] = sqrt((v[i + 1].first - v[i].first) * (v[i + 1].first - v[i].first)  + (v[i + 1].second - v[i].second) * (v[i + 1].second - v[i].second));
							lineSum[k] += lineTrack[k][str] * img[row][col];
							
			            }
			        }
				}
				if(lineSum[k] == 0) delta[k] = 0;
				else delta[k] = projection[j][k] / lineSum[k];
			}


			for(int k = 0; k < N; ++k) {
				for(unordered_map<string, double>::iterator itr = lineTrack[k].begin(); lineTrack[k].size() != 0 && itr != lineTrack[k].end(); ++itr) {
					string::size_type pos = itr->first.find(" ");
					string str = string(itr->first.begin(), itr->first.begin() + pos);
					int row = atoi(string(itr->first.begin(), itr->first.begin() + pos).c_str());
					int col = atoi(string(itr->first.begin() + pos + 1, itr->first.end()).c_str());
					
					pixelSum[row][col] += itr->second ;
					deltaMutilpyPixel[row][col] += itr->second * delta[k];
				}
				
			}
		}

		for(int i = 0; i < imgR; ++i) {
			for(int j = 0; j < imgC; ++j) {
				if(pixelSum[i][j] == 0) {
					img[i][j]  = 0;
				} else {
					double c = deltaMutilpyPixel[i][j] / pixelSum[i][j];
					img[i][j] *= c;
				}
			}
		}

		//save(i + 101, img);
	}

	int k = 0;
	for(int i = 0; i < imgWidth; ++i) {
		for(int j = 0; j < imgHeight; ++j) {
			output[k++] = img[j][i];
		}
	}
	
}
