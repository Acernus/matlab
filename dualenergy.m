function varargout = dualenergy(varargin)
% DUALENERGY MATLAB code for dualenergy.fig
%      DUALENERGY, by itself, creates a new DUALENERGY or raises the existing
%      singleton*.
%
%      H = DUALENERGY returns the handle to a new DUALENERGY or the handle to
%      the existing singleton*.
%
%      DUALENERGY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DUALENERGY.M with the given input arguments.
%
%      DUALENERGY('Property','Value',...) creates a new DUALENERGY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dualenergy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dualenergy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dualenergy

% Last Modified by GUIDE v2.5 14-Dec-2015 19:00:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dualenergy_OpeningFcn, ...
                   'gui_OutputFcn',  @dualenergy_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dualenergy is made visible.
function dualenergy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dualenergy (see VARARGIN)

% Choose default command line output for dualenergy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dualenergy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dualenergy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in readprojection.
function readprojection_Callback(hObject, eventdata, handles)
% hObject    handle to readprojection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
highprj = readbin('PROJ_MATCH_H2.BIN');
lowprj = readbin('PROJ_MATCH_L2.BIN');
sub = highprj - lowprj;
subplot(1, 3, 1); imshow(lowprj, [-0.02 4.43]);title('低能投影图');
subplot(1, 3, 2); imshow(highprj, [-0.01 3.4]);title('高能投影图');
subplot(1, 3, 3); imshow(sub, [0 1]);title('高低能差值图');


% --- Executes on button press in prjMatch.
function prjMatch_Callback(hObject, eventdata, handles)
% hObject    handle to prjMatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
highprj = readbin('PROJ_MATCH_H2.BIN');
lowprj = readbin('PROJ_MATCH_L2.BIN');
for i = 1 : 3
    highprj = [highprj(2:400,:);highprj(1,:)];
end
sub = highprj - lowprj;
subplot(1, 3, 1); imshow(lowprj, [-0.02 4.43]);title('低能投影图');
subplot(1, 3, 2); imshow(highprj, [-0.01 3.4]);title('高能投影图');
subplot(1, 3, 3); imshow(sub, [0 1]);title('高低能差值图');


% --- Executes on button press in prjMatching.
function prjMatching_Callback(hObject, eventdata, handles)
% hObject    handle to prjMatching (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;
plr = load('basic_lookuptable_b1.txt');
phr = load('basic_lookuptable_b2.txt');
prj1 = readbin('lookuptable_basic_decompositionprj1.bin');
prj2 = readbin('lookuptable_basic_decompositionprj2.bin');
lbd1 = readbin('lookuptable_basic_diff1.bin'); 
lbd2 = readbin('lookuptable_basic_diff2.bin'); 
lbd1 = reshape(lbd1, 1, 400 * 256);
lbd2 = reshape(lbd2, 1, 400 * 256);
subplot(2, 3, 1);
surf(B1, B2, plr);
axis tight;
colormap(hot);
shading interp;
title('低能查找表');
subplot(2, 3, 4);
surf(B1, B2, phr);
axis tight;
colormap(hot);
shading interp;
title('高能查找表');
subplot(2, 3, 2); imshow(prj1, []);title('C的投影分解图');
subplot(2, 3, 5); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 3, 3); plot(lbd1);title('高能误差分布图');
subplot(2, 3, 6); plot(lbd2);title('低能误差分布图');



% --- Executes on button press in energyDetail.
function energyDetail_Callback(hObject, eventdata, handles)
% hObject    handle to energyDetail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
low = load('newSpectrum80kV.dat');
lowx = zeros(1, 600);
highx = zeros(1, 600);
for i = 1 : 600
    lowx(i) = 1 + (i - 1) * 79 / 600;
    highx(i) = 1 + (i - 1) * 160 / 600;
end
high = load('newSpectrum160kV.dat');
subplot(2, 1, 1); plot(lowx, low); xlabel('keV');ylabel('Phatoms/mm^2');title('低能能谱');
subplot(2, 1, 2); plot(highx, high); xlabel('keV');ylabel('Phatoms/mm^2');title('高能能谱');

% --- Executes on button press in attenuation.
function attenuation_Callback(hObject, eventdata, handles)
% hObject    handle to attenuation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
low_al = load('0_80kv_Al.txt');
low_alx = low_al(:,1);
low_aly = low_al(:,2);
high_al = load('0_160kv_Al.txt');
high_alx = high_al(:,1);
high_aly = high_al(:,2);
low_c = load('0_80kv_Carbon.txt');
low_cx = low_c(:,1);
low_cy = low_c(:,2);
high_c = load('0_160kv_Carbon.txt');
high_cx = high_c(:,1);
high_cy = high_c(:,2);
subplot(2, 2, 1); plot(low_alx, low_aly);xlabel('keV');ylabel('cm^-1');title('Al的低能衰减系数');
subplot(2, 2, 2); plot(high_alx, high_aly);xlabel('keV');ylabel('cm^-1');title('Al的高能衰减系数');
subplot(2, 2, 3); plot(low_cx, low_cy);xlabel('keV');ylabel('cm^-1');title('C的低能衰减系数');
subplot(2, 2, 4); plot(high_cx, high_cy);xlabel('keV');ylabel('cm^-1');title('C的高能衰减系数');


% --- Executes on button press in backMatching.
function backMatching_Callback(hObject, eventdata, handles)
% hObject    handle to backMatching (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;
plr = load('lm_lookuptable_b1.txt');
phr = load('lm_lookuptable_b2.txt');
prj1 = readbin('lookuptable_lm_decompositionprj1.bin');
prj2 = readbin('lookuptable_lm_decompositionprj2.bin');
lbd1 = readbin('lookuptable_lm_diff1.bin'); 
lbd2 = readbin('lookuptable_lm_diff2.bin'); 
lbd1 = reshape(lbd1, 1, 400 * 256);
lbd2 = reshape(lbd2, 1, 400 * 256);
subplot(2, 3, 1);
surf(B1, B2, plr);
axis tight;
colormap(hot);
shading interp;
title('C的分解系数查找表');
subplot(2, 3, 4);
surf(B1, B2, phr);
axis tight;
colormap(hot);
shading interp;
title('Al的分解系数查找表');
subplot(2, 3, 2); imshow(prj1, []);title('C的投影分解图');
subplot(2, 3, 5); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 3, 3); plot(lbd1);title('高能误差分布图');
subplot(2, 3, 6); plot(lbd2);title('低能误差分布图');


% --- Executes on button press in lma.
function lma_Callback(hObject, eventdata, handles)
% hObject    handle to lma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

prj1 = readbin('lm_decomposition_prj1.bin');
prj2 = readbin('lm_decomposition_prj2.bin');
lbd1 = readbin('lm_diff1.bin'); 
lbd2 = readbin('lm_diff2.bin'); 
lbd1 = reshape(lbd1, 1, 400 * 256);
lbd2 = reshape(lbd2, 1, 400 * 256);
subplot(2, 2, 1); imshow(prj1, []);title('C的投影分解图');
subplot(2, 2, 3); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 2, 2); plot(lbd1);title('高能误差分布图');
subplot(2, 2, 4); plot(lbd2);title('低能误差分布图');


% --- Executes on button press in lmamatching.
function lmamatching_Callback(hObject, eventdata, handles)
% hObject    handle to lmamatching (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;
plr = load('lm_lookuptable_b1.txt');
phr = load('lm_lookuptable_b2.txt');
prj1 = readbin('combined_lookuptable_lm_decompositionprj1.bin');
prj2 = readbin('combined_lookuptable_lm_decompositionprj2.bin');
lbd1 = readbin('combined_lookuptable_lm_diff1.bin'); 
lbd2 = readbin('combined_lookuptable_lm_diff2.bin'); 
lbd1 = reshape(lbd1, 1, 400 * 256);
lbd2 = reshape(lbd2, 1, 400 * 256);
subplot(2, 3, 1);
surf(B1, B2, plr);
axis tight;
colormap(hot);
shading interp;
title('C的分解系数查找表');
subplot(2, 3, 4);
surf(B1, B2, phr);
axis tight;
colormap(hot);
shading interp;
title('Al的分解系数查找表');
subplot(2, 3, 2); imshow(prj1, []);title('C的投影分解图');
subplot(2, 3, 5); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 3, 3); plot(lbd1);title('高能误差分布图');
subplot(2, 3, 6); plot(lbd2);title('低能误差分布图');


% --- Executes on button press in fbpreconstruction.
function fbpreconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to fbpreconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

prj1 = readbin('lm_decomposition_prj1.bin');
prj2 = readbin('lm_decomposition_prj2.bin');
res1 = readbin('lm_FBP_1.bin');
res2 = readbin('lm_FBP_2.bin');
subplot(2, 2, 1); imshow(prj1, []);title('C的投影分解图');
subplot(2, 2, 2); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 2, 3); imshow(res1, []);title('重建的C的分解系数图');
subplot(2, 2, 4); imshow(res2, []);title('重建的Al的分解系数图');


% --- Executes on button press in mlem.
function mlem_Callback(hObject, eventdata, handles)
% hObject    handle to mlem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

prj1 = readbin('lm_decomposition_prj1.bin');
prj2 = readbin('lm_decomposition_prj2.bin');
res1 = readbin('lm_ml_em_img1.bin');
res2 = readbin('lm_ml_em_img2.bin');
subplot(2, 2, 1); imshow(prj1, []);title('C的投影分解图');
subplot(2, 2, 2); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 2, 3); imshow(res1, []);title('重建的C的分解系数图');
subplot(2, 2, 4); imshow(res2, []);title('重建的Al的分解系数图');


% --- Executes on button press in noisefbp.
function noisefbp_Callback(hObject, eventdata, handles)
% hObject    handle to noisefbp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
prj1 = imread('mark_fbp1.jpg');
prj2 = imread('mark_fbp2.jpg');
res1 = imread('fbp1.jpg');
res2 = imread('fbp2.jpg');

subplot(2, 2, 1); imshow(prj1, []);title('C的分解系数图');
subplot(2, 2, 3); imshow(prj2, []);title('Al的分解系数图');
subplot(2, 2, 2); imshow(res1, []);title('C的分解系数图噪声情况');
subplot(2, 2, 4); imshow(res2, []);title('Al的分解系数图噪声情况');


% --- Executes on button press in simulatemodel.
function simulatemodel_Callback(hObject, eventdata, handles)
% hObject    handle to simulatemodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
PicN = 512;
lowModelIntensity = [0.2 0.3 0.5 0.7 0.9];
highModelIntensity = [0.15 0.21 0.35 0.49 0.63];
lowModel = create(lowModelIntensity, PicN);
highModel = create(highModelIntensity, PicN);
subplot(1, 2, 1); imshow(lowModel, []);title('模拟低能衰减系数模型');
subplot(1, 2, 2); imshow(highModel, []);title('模拟高能衰减系数模型');

% --- Executes on button press in single.
function single_Callback(hObject, eventdata, handles)
% hObject    handle to single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
prj1 = readbin('single_decompositionprj1.bin');
prj2 = readbin('single_decompositionprj2.bin');
lbd1 = readbin('single_diff1.bin'); 
lbd2 = readbin('single_diff2.bin'); 
lbd1 = reshape(lbd1, 1, 400 * 256);
lbd2 = reshape(lbd2, 1, 400 * 256);
subplot(2, 2, 1); imshow(prj1, []);title('C的投影分解图');
subplot(2, 2, 3); imshow(prj2, []);title('Al的投影分解图');
subplot(2, 2, 2); plot(lbd1);title('高能误差分布图');
subplot(2, 2, 4); plot(lbd2);title('低能误差分布图');


% --- Executes on button press in simulitereconstruction.
function simulitereconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to simulitereconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);

R0 = 400; R1 = 0; M = 360; N = 600;
lowModelIntensity = [0.2 0.3 0.5 0.7 0.9];
highModelIntensity = [0.15 0.21 0.35 0.49 0.63];
%获取分解后的投影
g = getProjection(lowModelIntensity, R0, R1, M, N);
g1 = getProjection(highModelIntensity, R0, R1, M, N);
% In high voltage u1 = 1, u2 = 1.2, in low voltage u1 = 1.2, u2 = 1.5
lu1 = 0.2; lu2 = 0.3;
hu1 = 0.15; hu2 = 0.21;

[row, col] = size(g);
t1 = zeros(row, col);
t2 = zeros(row, col);
for i = 1 : row 
    for j = 1 : col
        t2(i, j) = (hu1 * g(i, j) - lu1 * g1(i, j)) / (hu1 * lu2 - lu1 * hu2);
        t1(i, j) = (g1(i, j) - t2(i, j) * hu2) / hu1;
    end
end
rt1 = imread('model_fbp1.jpg');
rt2 = imread('model_fbp2.jpg');
subplot(2, 3, 1); imshow(g, []);title('低能投影图');
subplot(2, 3, 4); imshow(g1, []);title('高能投影图');
subplot(2, 3, 2); imshow(t1, []);title('分解图');
subplot(2, 3, 5); imshow(t2, []);title('分解图');
subplot(2, 3, 3); imshow(rt1, []);title('重建图');
subplot(2, 3, 6); imshow(rt2, []);title('重建图');


% --- Executes on button press in unclefbp.
function unclefbp_Callback(hObject, eventdata, handles)
% hObject    handle to unclefbp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
prj1 = imread('mark_ml1.jpg');
prj2 = imread('mark_ml2.jpg');
res1 = imread('ml1.jpg');
res2 = imread('ml2.jpg');

subplot(2, 2, 1); imshow(prj1, []);title('C的分解系数图');
subplot(2, 2, 3); imshow(prj2, []);title('Al的分解系数图');
subplot(2, 2, 2); imshow(res1, []);title('C的分解系数图噪声情况');
subplot(2, 2, 4); imshow(res2, []);title('Al的分解系数图噪声情况');

% --- Executes on button press in lm_ml.
function lm_ml_Callback(hObject, eventdata, handles)
% hObject    handle to lm_ml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panhandle = uipanel(handles.uipanel);
axes('position',[2,0.714,129.6,29.929],'tag','imgdoc','box','off','NextPlot','replace','XGrid','off','YGrid','off','ZGrid','off', 'parent', panhandle);
prj1 = imread('mark_ml1.jpg');
prj2 = imread('mark_ml2.jpg');
res1 = imread('ml1.jpg');
res2 = imread('ml2.jpg');

subplot(2, 2, 1); imshow(prj1, []);title('C的分解系数图');
subplot(2, 2, 3); imshow(prj2, []);title('Al的分解系数图');
subplot(2, 2, 2); imshow(res1, []);title('C的分解系数图噪声情况');
subplot(2, 2, 4); imshow(res2, []);title('Al的分解系数图噪声情况');
