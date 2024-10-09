function varargout = GUIBEAM(varargin)
% GUIBEAM MATLAB code for GUIBEAM.fig
%      GUIBEAM, by itself, creates a new GUIBEAM or raises the existing
%      singleton*.
%
%      H = GUIBEAM returns the handle to a new GUIBEAM or the handle to
%      the existing singleton*.
%
%      GUIBEAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIBEAM.M with the given input arguments.
%
%      GUIBEAM('Property','Value',...) creates a new GUIBEAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIBEAM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIBEAM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIBEAM

% Last Modified by GUIDE v2.5 06-Jun-2014 19:44:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIBEAM_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIBEAM_OutputFcn, ...
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


% --- Executes just before GUIBEAM is made visible.
function GUIBEAM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIBEAM (see VARARGIN)

% Choose default command line output for GUIBEAM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIBEAM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIBEAM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global U B FFx FFy FBx FBy x_length y_length j x_element y_element e noo width a b x_node y_node element_number node_number d_x d_y X Y num A h g w kn D_matrix B_matrix u v Fx Fy nforce_number nBC_number fsurface_y fsurface_x fbody_x fbody_y
x_length=get(handles.edit1,'String');
x_length=str2num(x_length);
y_length=get(handles.edit2,'String');
y_length=str2num(y_length);
x_element=get(handles.edit3,'String');
x_element=str2num(x_element);
y_element=get(handles.edit4,'String');
y_element=str2num(y_element);
e=get(handles.edit5,'String');
e=str2num(e);
noo=get(handles.edit6,'String');
noo=str2num(noo);
width=get(handles.edit7,'String');
width=str2num(width);
element_number=2*x_element*y_element;
a=element_number;
x_node=1+x_element;
y_node=1+y_element;
node_number=x_node*y_node;
b=node_number;
d_x=x_length/x_element;
d_y=y_length/y_element;

 
X=zeros(node_number,1);
Y=zeros(node_number,1);
num=zeros(node_number,1);
A=zeros(a,1);
h=zeros(a,1);
g=zeros(a,1);
w=zeros(a,1);
kn=zeros(2*b);
h=zeros(element_number,1);
g=zeros(element_number,1);
w=zeros(element_number,1);
U=sym('U',[2*b 1]);
B=sym('B',[2*b 1]);
u=sym('u',[b 1]);
v=sym('v',[b 1]);
Fx=sym('Fx',[b 1]);
Fy=sym('Fy',[b 1]);
 j=[];
 nforce_number=[];
 nBC_number=[];
   fsurface_y=zeros(b,1);
    fsurface_x=zeros(b,1);
    fbody_x=zeros(b,1);
    fbody_y=zeros(b,1);
FFx=sym('FFx',[b 1]);
FFy=sym('FFy',[b 1]);
FBx=zeros(b,1);
FBy=zeros(b,1);
if get(handles.radiobutton1,'Value') == get(handles.radiobutton1,'Max')
     D_matrix=(e/((1+noo)*(1-noo)))*[1 noo 0;noo 1 0;0 0 (1-noo)/2];
 else
     D_matrix=(e/((1+noo)*(1-2*noo)))*[1-noo noo 0;noo 1-noo 0;0 0 (1-2*noo)/2];
 end
ii=1;
if get(handles.radiobutton5,'Value') == get(handles.radiobutton5,'Max')
for j=1:x_node
    for i=1:y_node
        num(ii,1)=(j-1)*y_node+i;
        X(num(ii))=(j-1)*d_x;
        Y(num(ii))=(i-1)*d_y;
        ii=ii+1;
        
    end
end
for i=1:element_number/2
    j1=floor(i/y_element);
    if(i/y_element)-round(i/y_element)==0
        j1=j1-1;
    end
    j2=ceil(i/y_element);
    h([2*i-1,2*i])=[i+j1;i+j1+y_node];
   g([2*i-1,2*i])=[i+(y_element+j2);i+y_element+j2+1];
   w([2*i-1,2*i])=[i+j2;i+j2];
  
end
else 
    for j=1:y_node
    for i=1:x_node
        num(ii,1)=(j-1)*x_node+i;
        X(num(ii))=(i-1)*d_x;
        Y(num(ii))=(j-1)*d_y;
        ii=ii+1;
        
    end
    end
%mohasebe ye node ha ye  dar bar girande ye har eleman 

for i=1:element_number/2
    j1=floor(i/x_element);
    if(i/x_element)-round(i/x_element)==0
        j1=j1-1;
    end
    j2=ceil(i/x_element);
    h([2*i-1,2*i])=[i+j1;i+j1];
   g([2*i-1,2*i])=[i+j2;i+(x_node+j2)];
   w([2*i-1,2*i])=[i+(x_node+j2);i+x_node+j1];
  
end
end
%%%%
 for i=1:element_number
      line([X(h(i,1),1) X(g(i,1),1)],[Y(h(i,1),1) Y(g(i,1),1) ],...
             'LineWidth',1,'Color',[0 0 1]);
    line([X(g(i,1),1) X(w(i,1),1)],[Y(g(i,1),1) Y(w(i,1),1) ],...
             'LineWidth',1,'Color',[0 0 1]);
    line([X(w(i,1),1) X(h(i,1),1)],[Y(w(i,1),1) Y(h(i,1),1) ],...
             'LineWidth',1,'Color',[0 0 1]);
   grid off
 end
 hold on
    
 for i=1:a
    A(i,1)=abs(0.5*det([1 X(h(i,1),1) Y(h(i,1),1);1 X(g(i,1),1) Y(g(i,1),1);1 X(w(i,1),1) Y(w(i,1),1)]));
     
     B_matrix{i}=(1/(2*A(i,1)))*[-(Y(w(i,1),1)-Y(g(i,1),1)) 0 -(Y(h(i,1),1)-Y(w(i,1),1)) 0 -(Y(g(i,1),1)-Y(h(i,1),1)) 0
         0 (X(w(i,1),1)-X(g(i,1),1)) 0 (X(h(i,1),1)-X(w(i,1),1)) 0 (X(g(i,1),1)-X(h(i,1),1))
         (X(w(i,1),1)-X(g(i,1),1)) -(Y(w(i,1),1)-Y(g(i,1),1)) (X(h(i,1),1)-X(w(i,1),1)) -(Y(h(i,1),1)-Y(w(i,1),1)) (X(g(i,1),1)-X(h(i,1),1)) -(Y(g(i,1),1)-Y(h(i,1),1))];
     k1{i}=(A(i)*width)*transpose(B_matrix{i})*D_matrix*B_matrix{i};
 end
for i=1:a
   
    K{i}=zeros(2*b);

    K{i}(2*h(i,1)-1:2*h(i,1),2*h(i,1)-1:2*h(i,1))=k1{i}(1:2,1:2);
    K{i}(2*h(i,1)-1:2*h(i,1),2*g(i,1)-1:2*g(i,1))=k1{i}(1:2,3:4);
    K{i}(2*h(i,1)-1:2*h(i,1),2*w(i,1)-1:2*w(i,1))=k1{i}(1:2,5:6);
    K{i}(2*g(i,1)-1:2*g(i,1),2*h(i,1)-1:2*h(i,1))=k1{i}(3:4,1:2);
    K{i}(2*g(i,1)-1:2*g(i,1),2*g(i,1)-1:2*g(i,1))=k1{i}(3:4,3:4);
    K{i}(2*g(i,1)-1:2*g(i,1),2*w(i,1)-1:2*w(i,1))=k1{i}(3:4,5:6);
    K{i}(2*w(i,1)-1:2*w(i,1),2*h(i,1)-1:2*h(i,1))=k1{i}(5:6,1:2);
    K{i}(2*w(i,1)-1:2*w(i,1),2*g(i,1)-1:2*g(i,1))=k1{i}(5:6,3:4);
    K{i}(2*w(i,1)-1:2*w(i,1),2*w(i,1)-1:2*w(i,1))=k1{i}(5:6,5:6);
end
for i=1:a
    kn=kn+K{i};
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1



% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2




function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  fsurface_x fsurface_y fbody_x fbody_y nBC_number j nforce_number b Fx Fy U B u v X Y XX YY h g w BB Bnew kn a B_matrix D_matrix sigma kornesh
if get(handles.radiobutton9,'Value')==1
   for i=1:b
    if X(i,1)==0
    nBC_number(i,1)=i;
    end
   end
   nBC_number=nonzeros(nBC_number);
   for  i=1:length(nBC_number)
       
        j(2*nBC_number(i,1)-1,1)=2*nBC_number(i,1)-1;
          u(nBC_number(i,1),1)=0;
      
    
     j(2*nBC_number(i,1),1)=2*nBC_number(i,1);
        v(nBC_number(i,1),1)=0; 
    
    if isreal(u(nBC_number(i,1),1))==1
        Fx(nBC_number(i,1),1)=Fx(nBC_number(i,1),1) ;

    end
   
     if isreal(v(nBC_number(i,1),1))==1
         Fy(nBC_number(i,1),1)=Fy(nBC_number(i,1),1);

     end
   end
   c=nonzeros(j);
else
    for i=1:length(nBC_number) 
j(2*nBC_number(i,1),1)=2*nBC_number(i,1);
       j(2*nBC_number(i,1)-1,1)=2*nBC_number(i,1)-1;
    end

c=nonzeros(j);
end
nBC_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node=zeros(b,1);
     for i=1:b
        node(i,1)=i;
     end
     set_1=setdiff(node,nBC_number);
     set_2=setdiff(set_1,nforce_number);
     
      for i=1:length(set_2)
          Fx(set_2(i,1),1)=0;
          Fy(set_2(i,1),1)=0;
      end
      if isempty(fsurface_x)==1
          for i=1:b
              fsurface_x(i,1)=0;
          end
      elseif  isempty(fsurface_y)==1
         for i=1:b
              fsurface_y(i,1)=0;
         end
      elseif isempty(fsurface_y)==1 && isempty(fsurface_y)==1
         for i=1:b
              fsurface_x(i,1)=0;
              fsurface_x(i,1)=0;
              
         end
         
      end
       if get(handles.checkbox7,'Value') == get(handles.checkbox7,'Min')
    for i=1:b
        fsurface_x(i,1)=0;
        fsurface_y(i,1)=0;
    end
        
    end
      
    for i=1:b
     B(2*i-1,1)=Fx(i,1)+fbody_x(i,1)+fsurface_x(i);
      
    B(2*i,1)=Fy(i,1)+fbody_y(i,1)+fsurface_y(i);
    
U(2*i-1,1)=u(i,1);
    
    U(2*i,1)=v(i,1);
    end
 B
   
cf=[];
for i=1:length(B)
    if isreal(B(i,1))==1
        cf(i,1)=i;
    end
end
cff=nonzeros(cf);
%tashkil e matris eslah shode ye niroo
Bnew=sym('B',[length(B) 1]);
z=[];
for i=1:length(cff)
    z(i,1)=0;
    for j=1:length(c)
    z(i,1)=z(i,1)-kn(cff(i,1),c(j,1))*U(c(j,1),1);
    end
end
for i=1:length(cff)
        Bnew(cff(i,1),1)=B(cff(i,1),1)+z(i,1);
        
end

%-----------------------------
% matrise kahesh yafte (k)
k=[];
k=kn;
for i=1:length(c)
        k(c(i,1),:)=0;
        k(:,c(i,1))=0;
       k(c(i,1),c(i,1))=1;
end
 %amaliyat baraye mohasebe ye jabe jaE ha
     nodenum=zeros(2*b,1);
     for i=1:2*b
        nodenum(i,1)=i;
     end
    set=setdiff(nodenum,c);
    UU=sym('UU',[length(set) 1]);
    BB=sym('BB',[length(cff) 1]);
     for i=1:length(cff)
        BB(i,1)=Bnew(cff(i,1),1);
    end
   UU=inv(k(cff,cff))*BB;
   for i=1:length(UU)
       U(set(i,1),1)=UU(i,1);
   end
   Bnew=kn*U;
  
    
    %-----------------------------------------------------
   %amaliate mohasbe niroo
    U=sym2poly(U)
     Bnew=sym2poly(Bnew)
   % axis([-1.1*x_length 1.1*x_length -1.1*y_length 1.1*y_length])
        for i=1:a
        XX(h(i,1),1)=X(h(i,1),1)+U(2*h(i,1)-1,1);
        YY(h(i,1),1)=Y(h(i,1),1)+U(2*h(i,1),1);
        XX(g(i,1),1)=X(g(i,1),1)+U(2*g(i,1)-1,1);
        YY(g(i,1),1)=Y(g(i,1),1)+U(2*g(i,1),1);
        XX(w(i,1),1)=X(w(i,1),1)+U(2*w(i,1)-1,1);
        YY(w(i,1),1)=Y(w(i,1),1)+U(2*w(i,1),1);
   line([XX(h(i,1),1) XX(g(i,1),1)],[YY(h(i,1),1) YY(g(i,1),1) ],...
             'LineWidth',0.5,'Color',[1 0 0]);
    line([XX(g(i,1),1) XX(w(i,1),1)],[YY(g(i,1),1) YY(w(i,1),1) ],...
             'LineWidth',0.5,'Color',[1 0 0]);
    line([XX(w(i,1),1) XX(h(i,1),1)],[YY(w(i,1),1) YY(h(i,1),1) ],...
             'LineWidth',0.5,'Color',[1 0 0]);
        end
        sigma=[];
        kornesh=[];
    for i=1:a
        sigma{i}=D_matrix*B_matrix{i}*U([2*h(i,1)-1 2*h(i,1) 2*g(i,1)-1 2*g(i,1) 2*w(i,1)-1 2*w(i,1)],1);
        kornesh{i}=B_matrix{i}*U([2*h(i,1)-1 2*h(i,1) 2*g(i,1)-1 2*g(i,1) 2*w(i,1)-1 2*w(i,1)],1);
    end
   
 
function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
global X Y
for  i=1:length(X)
        text(X(i) , Y(i),num2str(i))
 end

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
     %namgozari bar rooye eleman ha be tori ke dar vasate elemane cst
    %shomarande gharar migirad
    global X Y h g w element_number
    for i=1:element_number
       x_ave=(X(h(i))+X(g(i))+X(w(i)))/3;
     
          y_ave=(Y(h(i))+Y(g(i))+Y(w(i)))/3;
        text(x_ave,y_ave,num2str(i))
    end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   v u nBC_number   
z=get(handles.edit8,'string');
z=str2double(z);

for i=1:z
    if get(handles.popupmenu2,'Value')==i
        b1=get(handles.edit12,'string');
        b1=str2double(b1);
       nBC_number(i,1)=b1;
       u(b1)=get(handles.edit10,'String');
     v(b1)=get(handles.edit11,'String');
       
       
    end
end

nBC_number

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit1,'string',' ')
set(handles.edit2,'string',' ')
set(handles.edit3,'string',' ')
set(handles.edit4,'string',' ')
set(handles.edit5,'string',' ')
set(handles.edit6,'string',' ')
set(handles.edit7,'string',' ')
set(handles.edit8,'string',' ')
set(handles.edit12,'string',' ')
set(handles.edit10,'string',' ')
set(handles.edit11,'string',' ')
set(handles.edit21,'string',' ')
set(handles.edit18,'string',' ')
set(handles.edit20,'string',' ')
set(handles.edit19,'string',' ')
hold off
plot(handles.axes1,[0],[0])

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z=get(handles.edit8,'String');
z=str2double(z);
a=[1];
for i=2:z
a=[a;i];
end
set(handles.popupmenu2,'String',a)


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global  nforce_number Fx Fy
z=get(handles.edit21,'string');
z=str2double(z);
for i=1:z
    if get(handles.popupmenu4,'Value')==i
        b2=get(handles.edit18,'string');
        b2=str2double(b2);
       Fx(b2)=get(handles.edit20,'string');
       Fy(b2)=get(handles.edit19,'string');
       nforce_number(b2)=b2;
    end
end 

nforce_number=nonzeros(nforce_number)
Fx
Fy
    
% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z=get(handles.edit21,'String');
z=str2double(z);
a1=[1];
for i=2:z
a1=[a1;i];
end
set(handles.popupmenu4,'String',a1)



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
global a b fsurface_x fsurface_y

if get(handles.checkbox7,'Value') == get(handles.checkbox7,'Max')
    set(handles.text33,'Enable','on')
    set(handles.text32,'Enable','on')
    set(handles.text34,'Enable','on')
    set(handles.edit24,'Enable','on')
    set(handles.edit23,'Enable','on')
    set(handles.edit25,'Enable','on')
    set(handles.text31,'Enable','on')
    set(handles.edit22,'Enable','on')
    set(handles.popupmenu5,'Enable','on')
    set(handles.pushbutton11,'Enable','on')
    set(handles.pushbutton10,'Enable','on')
else
     set(handles.text33,'Enable','off')
    set(handles.text32,'Enable','off')
    set(handles.text34,'Enable','off')
    set(handles.edit24,'Enable','off')
    set(handles.edit23,'Enable','off')
    set(handles.edit25,'Enable','off')
    set(handles.text31,'Enable','off')
    set(handles.edit22,'Enable','off')
    set(handles.popupmenu5,'Enable','off')
    set(handles.pushbutton11,'Enable','off')
    set(handles.pushbutton10,'Enable','off')
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global surfnum Y X fsurface_x fsurface_y
syms x
for i=1:surfnum
    if get(handles.popupmenu5,'Value')==i
             surfnode_first=get(handles.edit23,'String');
             surfnode_first=str2double(surfnode_first);
           surfnode_sec=get(handles.edit25,'String');
           surfnode_sec=str2double(surfnode_sec);
           l(i)=sqrt((X(surfnode_first)-X(surfnode_sec))^2+(Y(surfnode_first)-Y(surfnode_sec))^2);
            
           wx=get(handles.edit22,'String');
         
           Rx(i)=int(wx,x,0,l(i));
           xbar(i,1)=int(sym(x*wx)/Rx,x,0,l(i));
           P1_surface(i)=Rx(i)*(l(i)-xbar(i ))/l(i);
       P2_surface(i)=Rx(i)*xbar(i,1)/l(i);
       if Y(surfnode_first)==Y(surfnode_sec)
          
            fsurface_y(surfnode_first)=P1_surface(i);
      
           fsurface_y(surfnode_sec)=P2_surface(i);
       elseif X(surfnode_first)==X(surfnode_sec)
           
            fsurface_x(surfnode_first)=P1_surface(i);
      
           fsurface_x(surfnode_sec)=P2_surface(i);
       end
    end
end
fsurface_x
fsurface_y
% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global surfnum
surfnum=get(handles.edit24,'String');
surfnum=str2double(surfnum);
a=[1];
for i=2:surfnum
a=[a;i];
end

set(handles.popupmenu5,'string',a)


function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8

if get(handles.checkbox8,'Value') == get(handles.checkbox8,'Max')
    set(handles.text35,'Enable','on')
    set(handles.text36,'Enable','on')
    set(handles.pushbutton12,'Enable','on')
    set(handles.edit26,'Enable','on')
    set(handles.edit27,'Enable','on')
else
        set(handles.text35,'Enable','off')
    set(handles.text36,'Enable','off')
    set(handles.pushbutton12,'Enable','off')
    set(handles.edit26,'Enable','off')
    set(handles.edit27,'Enable','off')
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global A width a  w g h body_x body_y fbody_x fbody_y
body_x=get(handles.edit26,'string');
body_x=str2double(body_x);
body_y=get(handles.edit27,'string');
body_y=str2double(body_y);
for i=1:a
    fbody_x(h(i,1),1)= fbody_x(h(i,1),1)+((A(i,1)*width)/3)*body_x;
          fbody_x(g(i,1),1)= fbody_x(g(i,1),1)+((A(i,1)*width)/3)*body_x;
          fbody_x(w(i,1),1)=fbody_x(w(i,1),1)+((A(i,1)*width)/3)*body_x;
          fbody_y(h(i,1),1)=fbody_y(h(i,1),1)+ ((A(i,1)*width)/3)*body_y;
          fbody_y(g(i,1),1)= fbody_y(g(i,1),1)+((A(i,1)*width)/3)*body_y;
          fbody_y(w(i,1),1)= fbody_y(w(i,1),1)+((A(i,1)*width)/3)*body_y;
end
fbody_x
fbody_y


% --- Executes during object creation, after setting all properties.
function checkbox7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global   a b  Bnew U sigma kornesh
set(handles.uipanel10,'Visible','on')
set(handles.axes1,'visible','off')
u=[];
v=[];
Fx=[];
Fy=[];
for i=1:b
   
    Fx(i)=Bnew(2*i-1);
    Fy(i)=Bnew(2*i);
    u(i)=U(2*i-1);
    v(i)=U(2*i);
    
end
for i=1:a
    sigma_x(i,1)=sigma{i}(1);
    sigma_y(i,1)=sigma{i}(2);
    taw_xy(i,1)=sigma{i}(3);
    epsilon_x(i,1)=kornesh{i}(1);
    epsilon_y(i,1)=kornesh{i}(2);
    gama_xy(i,1)=kornesh{i}(3);
end

if   get(handles.popupmenu6,'Value')==1
set(handles.listbox1,'string',u)
elseif get(handles.popupmenu6,'Value')==2
    set(handles.listbox1,'string',v)
elseif get(handles.popupmenu6,'Value')==3
set(handles.listbox1,'string',Fx)
elseif get(handles.popupmenu6,'Value')==4
    
set(handles.listbox1,'string',Fy)
elseif get(handles.popupmenu6,'Value')==5
    
set(handles.listbox1,'string',sigma_x)

elseif get(handles.popupmenu6,'Value')==6
    
set(handles.listbox1,'string',sigma_y)

elseif get(handles.popupmenu6,'Value')==7
    
set(handles.listbox1,'string',taw_xy)

elseif get(handles.popupmenu6,'Value')==8
    
set(handles.listbox1,'string',epsilon_x)

elseif get(handles.popupmenu6,'Value')==9
    
set(handles.listbox1,'string',epsilon_y)

elseif get(handles.popupmenu6,'Value')==10
    
set(handles.listbox1,'string',gama_xy)
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel10,'Visible','off')
set(handles.axes1,'visible','on')


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9
set(handles.uipanel2,'Visible','off')

% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10
set(handles.uipanel2,'Visible','on')
