function varargout = two_layer_immersed_v005(varargin)
% Shull research group. Original author: Chyi-Huey Joshua Yeh
%This version is currently in beta mode
% TWO_LAYER_IMMERSED_V005 MATLAB code for two_layer_immersed_v005.fig
%      TWO_LAYER_IMMERSED_V005, by itself, creates a new TWO_LAYER_IMMERSED_V005 or raises the existing
%      singleton*.
%
%      H = TWO_LAYER_IMMERSED_V005 returns the handle to a new TWO_LAYER_IMMERSED_V005 or the handle to
%      the existing singleton*.
%
%      TWO_LAYER_IMMERSED_V005('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWO_LAYER_IMMERSED_V005.M with the given input arguments.
%
%      TWO_LAYER_IMMERSED_V005('Property','Value',...) creates a new TWO_LAYER_IMMERSED_V005 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before two_layer_immersed_v005_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to two_layer_immersed_v005_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help two_layer_immersed_v005

% Last Modified by GUIDE v2.5 04-Aug-2017 16:21:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @two_layer_immersed_v005_OpeningFcn, ...
                   'gui_OutputFcn',  @two_layer_immersed_v005_OutputFcn, ...
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


% --- Executes just before two_layer_immersed_v005 is made visible.
function two_layer_immersed_v005_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to two_layer_immersed_v005 (see VARARGIN)

% Choose default command line output for two_layer_immersed_v005
handles.output = hObject;

%set default values or settings
handles=home_state(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes two_layer_immersed_v005 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles=home_state(handles)
%This function contains all of the default settings of the GUI.
cla(handles.axes1);cla(handles.axes2);
filepath=[];
if isfield(handles,'din')==true
    if isfield(handles.din,'filepath')==true
        filepath=handles.din.filepath;
    end
end
handles.din=[];
if isempty(filepath)==0
    handles.din.filepath=filepath;
end
%set default values or settings
warning('off','MATLAB:legend:IgnoringExtraEntries');%suppress warning messages associated with havein extra legend entries
handles.din.marker_color=[{[0 0 1]} {[1 0 0]} {[0 0.5 0]} {[0.5 0.5 0]} {[0 0.5 0.5]} {[0.5 0 0.5]}];%set marker colors
handles.din.marker_style=[{['+']} {['o']} {['s']} {['^']} {['v']} {['x']}];%set marker style
if isfield(handles.din,'filepath')==false
    handles.din.filepath=pwd;%set default working path directory
end
handles.din.constants.f1=5e6;%fundamental resonance frequency in Hz
handles.din.constants.zq=8.84e6; %load impedance of quartz, kg/m^2-s
handles.din.rawdata.immersed_freq=[];
handles.din.rawdata.immersed_diss=[];
handles.din.simulate=0;
handles.din.bare_flag=0;%this is flag, showing that the bare cyrstal data was not loaded
handles.din.bare_path=pwd;%create a fieldname for the path destination of the bare crystal file
handles.din.n1=3;%harmonic order in which solution will be calculated at
handles.din.fit_factor_range=6;%this factor influences the width or freq range in which the raw spectras will be fitted
handles.prefs.num_peaks=[1,0,1,0,1,0,1,0,1,0,1].*3;%set default option for the number of peaks the program will maximally fit
handles.prefs.sensitivity(1)=5e-3;%peak sensitivity factor for the 1st harmonic
handles.prefs.sensitivity(3)=5e-3;%peak sensitivity factor for the 3rd harmonic
handles.prefs.sensitivity(5)=5e-3;%peak sensitivity factor for the 5th harmonic
handles.prefs.sensitivity(7)=5e-3;%peak sensitivity factor for the 7th harmonic
handles.prefs.sensitivity(9)=5e-3;%peak sensitivity factor for the 9th harmonic
handles.prefs.sensitivity(11)=5e-3;%peak sensitivity factor for the 11th harmonic
handles.prefs.peak_min(1)=.2;%min. peak finding threshold for the 1st harmonic
handles.prefs.peak_min(3)=.2;%min. peak finding threshold for the 3rd harmonic
handles.prefs.peak_min(5)=.2;%min. peak finding threshold for the 5th harmonic
handles.prefs.peak_min(7)=.2;%min. peak finding threshold for the 7th harmonic
handles.prefs.peak_min(9)=.2;%min. peak finding threshold for the 9th harmonic
handles.prefs.peak_min(11)=.2;%min. peak finding threshold for the 11th harmonic
handles.prefs.simul_peak=1;
data=zeros(6,6);
for dum=1:12
    for dum1=1:12
        if dum1==dum&&mod(dum,2)==1
            data(dum,dum1)=(dum*2-1)*10;
        elseif dum1==dum&&mod(dum,2)==0
            data(dum,dum1)=((dum-1)*2-1)*10;
        end
    end
end
set(handles.error_values,'userdata',data);
handles.din.contour.label=0;%default flag for the harm and diss ratios
handles.din.guess_label=1;%default flag for the guess parameters
handles.din.guess_label2=1;%default flag for the guess parameters
handles.din.qcm_label=0;%default plaf for the qcm map
handles.din.contour.res=150;%define resolution for the contour plots
handles.din.lb=[0 0 0];%set lowerbound guess values (d2lam, phi, drho)
handles.din.ub=[1 90 10];%set upperbound guess values (d2lam, phi, drho)
handles.din.contour.harm_ratio_calc=ones(handles.din.contour.res,handles.din.contour.res);%allocate harm_ratio_calc
handles.din.contour.diss_ratio_calc=ones(handles.din.contour.res,handles.din.contour.res);%allocate diss_ratio_calc
handles.din.qcm_map=ones(handles.din.contour.res,handles.din.contour.res);%allocate the matrix for the qcm map
handles.din.grho_state=1;
handles.din.d2lam_state=0;
handles.din.fit_options=optimset('Display','off','tolfun',1e-11,'tolx',1e-11,'Maxiter',1e8,'MaxFunEvals',1e8);
set(handles.harm1,'value',1);
set(handles.harm3,'value',1);
set(handles.harm5,'value',1);
set(handles.radio5,'value',1);
set(handles.dgliq,'string',1000);%set default value for dgliq
set(handles.dgliq_harm,'string',3);%set default value for dgliq_harm
set(handles.solve_inc,'string',1);%set default value for the increment associated with the solve all functionality
set(handles.auto_solve,'value',0);
set(handles.edit_drho,'value',0);
set(handles.filename_txt,'string','<filename>','tooltipstring','');
set(handles.bare_table,'data',[]);
set(handles.uitable2,'data',[]);
set(handles.uitable3,'data',[]);
set(handles.status,'string','Status: Ready...','backgroundcolor','k','foregroundcolor','r');
set(handles.initial_time,'string',1);
set(handles.bare_xtal,'tooltipstring','');
set(handles.figure1,'units','pixels');
dcm=datacursormode(handles.figure1);
set(dcm,'UpdateFcn',{@myupdatefcn,handles},'enable','off');
for dum=1:1:6
    name=['harmchoice',num2str(dum)];
    set(handles.(name),'string',[{'1,3:1'},{'1,3:3'},{'1,5:1'},{'1,5:5'},{'3,5:3'},{'3,5:5'},...
        {'1,7:7'},{'3,7:7'},{'5,7:7'},{'1,7:1'},{'3,7:3'},{'5,7:5'},{'1:5,3'},{'5,3:3'},{'3,9:3'},...
        {'5,9:5'},{'5,9:9'},{'3,11:11'},{'3,11:3'},{'7,9:7'},{'7,9:9'},{'9,11:9'},{'9,11:11'},...
        {'7,11:7'},{'7,11:11'},{'3,5:7'},{'1,3:5'},{'3,5:1'}]);
end%for dum=1:1:6
set(handles.contour_choice,'string',get(handles.(name),'string'));
set(handles.bare_table,'ColumnName',{'<HTML><b>f (Hz)</b></html>','<HTML><b> &Gamma (Hz)</b></html>'});
set(handles.uitable2,'ColumnName',{'<HTML><b>d&rho (<sup>g</sup>&frasl;<sub>m</sub><font size="2">2</font>)</b></html>',...
    '<HTML><b>&#8739;G*&#8739;&rho (<sup>Pa-g</sup>&frasl;<sub>cm</sub><font size="2">3</font>)</b></html>',...
    '<HTML><b>&phi (deg.)</b></html>','<HTML><b><sup>d</sup>&frasl;<sub>&lambda;n</sub></b></html>',...
    '<HTML><b>&lambda;&rho; (<sup>&mu;m-g</sup>&frasl;<sub>cm</sub><font size="2">3</font>)</b></html>','<HTML><b><sup>&Delta;f</sup>&frasl;<sub>&Delta;fsn</sub></b></html>',...
    '<HTML><b><sup>&Delta;&Gamma;</sup>&frasl;<sub>&Delta;fsn</sub></b></html>','<HTML><b>cf</b></html>',...
    '<HTML><b><sup>&delta;</sup>/<sub>d</sub></b></html>'});
set(handles.uitable3,'ColumnName',{'<HTML><b>raw &Delta;f (Hz)</b></html>',...
    '<HTML><b>pred. &Delta;f (Hz)</b></html>','<HTML><b>raw &Delta;&Gamma; (Hz)</b></html>',...
    '<HTML><b>pred. &Delta;&Gamma; (Hz)</b></html>'});
set(handles.uitable4,'ColumnName',{'variables','guess values'},...
    'data',[{'<HTML><b>d&rho; (g/m<sup><font size=2>2</font></sup>)</b></html>'},{1};...
    {'<HTML><b>&phi; (deg.)</b></html>'},{45};...
    {'<HTML><body text=#FF0000><b>&#8739;G*&#8739;&rho; (Pa-g/cm<sup><font size=2>3</font></sup>)</body></b></html>'},{1e6};...
    {'<HTML><b>d/&lambda;</b></html>'},{0.02}]);
tt=text(0.25,0.25,'$\bf\Delta\Gamma$','parent',handles.dumaxes,'verticalalignment','bottom','fontweight','bold');
set(tt,'interpreter','latex');
set(handles.dumaxes,'xtick',[],'ytick',[],'box','off','xcolor','w','ycolor','w');
temp.Indices=[1,2];
uitable4_CellEditCallback(handles.uitable4, temp, handles);


% --- Outputs from this function are returned to the command line.
function varargout = two_layer_immersed_v005_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACK FUNCTION IN THE TOP "TOOLBAR" PANEL
function load_Callback(hObject, eventdata, handles)
%This function loads the frequency shift data (shift in the location of the
%resonance peak and the shift in the HMHW of the resonance peak).
handles=home_state(handles);%set the GUI state back to its default state
guidata(handles.figure1,handles);
set(handles.status,'string','Status: GUI handles have been reset to its default state!','backgroundcolor','k','foregroundcolor','r');
disp('GUI handles have been reset to its default state!');
handles.din.harmtot=active_harm(handles);%get the active harmonics
try%prompt the user to import the frequency shift values    
    [filename,filepath,~]=uigetfile('.mat','Load frequency shift data',handles.din.filepath);%get the filepath and filename
    if isstr(filename)==0
        return
    end
    set(handles.status,'string','Status: Importing...','backgroundcolor','k','foregroundcolor','r');
    handles=import_data(handles,filename,filepath);%import the data
    handles.din.filepath=filepath;
    handles.din.bare_path=filepath;
    handles.raw=[];
    set(hObject,'tooltipstring',['<html>Data loaded on ',datestr(clock),'<br />Filename: ',filename,'<br />Filepath: ',filepath,'</html>']);
    disp('Data successfully imported!');
    disp([filepath,filename]);
    set(handles.status,'string',['Status: Data successfully imported (',[filepath(1:10),'...',filename],')!'],'backgroundcolor','k','foregroundcolor','r');
    guidata(hObject, handles);
catch err_msg
    assignin('base','err_msg',err_msg);
    disp('Import failed!');
    set(handles.status,'string','Status: Import failed!','backgroundcolor','r','foregroundcolor','k');
end%try

function handles=import_data(handles,filename,filepath)
%this function loads the data and store it in the handles structure
handles.din.rawdata=load([filepath,filename]);%load the dataset
set(handles.filename_txt,'string',filename,'tooltipstring',[filepath,filename]);
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
handles=bare_crystal(handles);%output the reference crystal readings
timepoints=handles.din.rawdata.freq_shift(:,1);%xtract out the timepoints
ind=find(isnan(timepoints)==0);
timepoints=timepoints(ind);
initial_time_Callback(handles.initial_time, 1, handles);
set(handles.extent,'string',length(timepoints));%set end point to last timepoint
extent_Callback(handles.extent, 1, handles);
%update the handles structure
handles.din.filename=filename;
handles.din.filepath=filepath;   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACK FUNCTION IN THE PLOTTING OPTIONS PANEL

function plot1_options_Callback(~, ~, handles)
%This function allows the user to access the axes options for
%handles.axes1. See MATLAB function instpect for details.
inspect(handles.axes1);
function plot2_options_Callback(~, ~, handles)
%This function allows the user to access the axes options for
%handles.axes2. see MATLAB function for details.
inspect(handles.axes2);
function time_units_Callback(hObject, ~, handles)
%This function allows the user to adjust the units for the time axis
handles=plot_fcn(handles,handles.din.rawdata);
initial_time_Callback(handles.initial_time, 1, handles);
extent_Callback(handles.extent, 1, handles);
guidata(hObject,handles);
function norm_time_Callback(hObject, ~, handles)
%this function allows the user to choose whether or not to set the first
%timepoint as 0
handles=plot_fcn(handles,handles.din.rawdata);
initial_time_Callback(handles.initial_time, 1, handles);
extent_Callback(handles.extent, 1, handles);
guidata(hObject,handles);
function linkx_Callback(hObject, ~, handles)
%this function allows the user to choose whether or not to link the x-axis
%together.
handles=plot_fcn(handles,handles.din.rawdata);
initial_time_Callback(handles.initial_time, 1, handles);
extent_Callback(handles.extent, 1, handles);
guidata(hObject,handles);

% --- Executes on button press in MS.
function MS_Callback(hObject, eventdata, handles)
clf(figure(99));
set(figure(99),'menubar','none','name','Gen. Mat. Space Plots',...
    'numbertitle','off','color','w');
f=figure(99);
f.Position(3)=300;
handles2.drho_min=uicontrol('style','edit','units','normalized','tooltipstring','g/m^2',...
    'position',[0.35 0.901 0.3 0.05],'string',0.1);
handles2.drho_max=uicontrol('style','edit','units','normalized','tooltipstring','g/m^2',...
    'position',[0.35 0.851 0.3 0.05],'string',5);
handles2.phi_min=uicontrol('style','edit','units','normalized','tooltipstring','deg.',...
    'position',[0.35 0.751 0.3 0.05],'string',0);
handles2.phi_max=uicontrol('style','edit','units','normalized','tooltipstring','deg.',...
    'position',[0.35 0.701 0.3 0.05],'string',20);
handles2.grho_min=uicontrol('style','edit','units','normalized','tooltipstring','Pa-g/cm^3',...
    'position',[0.35 0.601 0.3 0.05],'foregroundcolor','r');
handles2.grho_max=uicontrol('style','edit','units','normalized','tooltipstring','Pa-g/cm^3',...
    'position',[0.35 0.551 0.3 0.05],'foregroundcolor','r');
handles2.d2lam_min=uicontrol('style','edit','units','normalized','tooltipstring','normlized thickness',...
    'position',[0.35 0.451 0.3 0.05],'string',0.02);
handles2.d2lam_max=uicontrol('style','edit','units','normalized','tooltipstring','normlized thickness',...
    'position',[0.35 0.401 0.3 0.05],'string',0.5);
handles2.delf_max=uicontrol('style','edit','units','normalized','tooltipstring','Hz',...
    'position',[0.35 0.301 0.3 0.05],'string',1e5);
handles2.delg_max=uicontrol('style','edit','units','normalized','tooltipstring','Hz',...
    'position',[0.35 0.251 0.3 0.05],'string',2e4);
handles2.drho_min_txt=uicontrol('style','text','units','normalized','string','drho min:',...
    'position',[0.05 0.9 0.3 0.05]);
handles2.drho_max_txt=uicontrol('style','text','units','normalized','string','drho max:',...
    'position',[0.05 0.85 0.3 0.05]);
handles2.phi_min_txt=uicontrol('style','text','units','normalized','string','phi min:',...
    'position',[0.05 0.75 0.3 0.05]);
handles2.phi_max_txt=uicontrol('style','text','units','normalized','string','phi max:',...
    'position',[0.05 0.7 0.3 0.05]);
handles2.grho_min_txt=uicontrol('style','text','units','normalized','string','grho min:',...
    'position',[0.05 0.6 0.3 0.05],'foregroundcolor','r');
handles2.grho_max_txt=uicontrol('style','text','units','normalized','string','grho max:',...
    'position',[0.05 0.55 0.3 0.05],'foregroundcolor','r');
handles2.d2lam_min_txt=uicontrol('style','text','units','normalized','string','d2lam1 min:',...
    'position',[0.05 0.45 0.3 0.05]);
handles2.d2lam_max_txt=uicontrol('style','text','units','normalized','string','d2lam1 max:',...
    'position',[0.05 0.4 0.3 0.05]);
handles2.delf_max_txt=uicontrol('style','text','units','normalized','string','\Deltaf max:',...
    'position',[0.05 0.3 0.3 0.05]);
handles2.delg_max_txt=uicontrol('style','text','units','normalized','string','\Deltag max:',...
    'position',[0.05 0.25 0.3 0.05]);
handles2.plot_button=uicontrol('style','pushbutton','units','normalized','string','Plot',...
    'fontsize',10,'fontweight','bold','backgroundcolor',[0.7 0.7 0.7],'position',[0.025 0.005 0.2 0.05]);
handles2.uibuttongrp=uibuttongroup('units','normalized','position',[0.025 0.075 0.4 0.15],'backgroundcolor','w');
handles2.one_layer=uicontrol('style','radiobutton','units','normalized','string','1-layer',...
    'position',[0.05 0.5 0.8 0.45],'backgroundcolor','w','parent',handles2.uibuttongrp);
handles2.two_layer=uicontrol('style','radiobutton','units','normalized','string','2-layer',...
    'tooltipstring','<html>impedance calc. of the semiinfinite layer is based on the dissipation values <br/>listed in the "Calculations options" panel in the main GUI',...
    'position',[0.05 0 0.8 0.45],'backgroundcolor','w','parent',handles2.uibuttongrp);
set(findall(f,'style','text'),'backgroundcolor','w','horizontalalignment','right',...
    'fontweight','bold','fontsize',10);
set(handles2.grho_min,'callback',{@d2lam_calc,handles,handles2});
set(handles2.grho_max,'callback',{@d2lam_calc,handles,handles2});
set(handles2.d2lam_min,'callback',{@grho_calc2,handles,handles2});
set(handles2.d2lam_max,'callback',{@grho_calc2,handles,handles2});
set(handles2.drho_min,'callback',{@MS_recalc,handles,handles2});
set(handles2.drho_max,'callback',{@MS_recalc,handles,handles2});
set(handles2.phi_min,'callback',{@MS_recalc,handles,handles2});
set(handles2.phi_max,'callback',{@MS_recalc,handles,handles2});
set(handles2.plot_button,'callback',{@MS_plot,handles,handles2});
grho_calc2(1,1,handles,handles2);
guidata(handles2.plot_button,handles2);

function MS_plot(~,~,handles,handles2)
handles=guidata(handles.figure1);%refresh the handles structure
drho=[str2double(handles2.drho_min.String) str2double(handles2.drho_max.String)];
phi=[str2double(handles2.phi_min.String) str2double(handles2.phi_max.String)];
grho=[str2double(handles2.grho_min.String) str2double(handles2.grho_max.String)];
d2lam=[str2double(handles2.d2lam_min.String) str2double(handles2.d2lam_max.String)];
del_fg=[str2double(handles2.delf_max.String) str2double(handles2.delg_max.String)];
f1=handles.din.constants.f1;%fundamentalr resonance freq in Hz
zq=handles.din.constants.zq;% the quartz load impedance
%get the harmonic combinations in which the calcualtions will be based off
%of
harm_value=handles.contour_choice.Value;
harm_string=handles.contour_choice.String;
MS_harm=harm_string{harm_value};
index1=find(MS_harm==',');
index2=find(MS_harm==':');
n1=str2double(MS_harm(1:index1-1));
n2=str2double(MS_harm(index1+1:index2-1));
n3=str2double(MS_harm(index2+1:end));
%run the calcualtions
%calculate the liquid load impedance
if handles2.two_layer.Value==1
    dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
    dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance at dgliq_harm
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at
else
    dgliq_n=0;%diss shift used for liq load impedance calc
    dfliq_n=0;%freq shift used for liq load impedance calc
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance at dgliq_harm
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at
end
unique_n=unique([n1 n2 n3]);
for dum=1:length(unique_n)
    name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
    handles.din.solved.(name).drho=drho;
    handles3=calc_viscoelastic_par(handles,unique_n(dum),f1,z_star_liq_n,zq,[d2lam phi drho],dgliq_harm,name,n1);%calculate the viscoelastic parameters
end
clf(figure(100));
f=figure(100);
set(f,'numbertitle','off');
keyboard


function MS_recalc(~,~,handles,handles2)
temp=findall(figure(99),'style','text','foregroundcolor','r');
if strcmp(temp(1).String,'grho min:')==1||strcmp(temp(1).String,'grho max:')==1
    grho_calc2(0,0,handles,handles2);
elseif strcmp(temp(1).String,'d2lam1 min:')==1||strcmp(temp(1).String,'d2lam1 max:')==1
    d2lam_calc(0,0,handles,handles2);
end

function grho_calc2(~,~,handles,handles2)
%minimum calculation
handles2.phi_val_min=str2double(handles2.phi_min.String);%phi (deg.)
handles2.drho_val_min=str2double(handles2.drho_min.String)./1000;%areal mass (kg/m^2)
handles2.d2lam_val_max=str2double(handles2.d2lam_max.String);
f1=handles.din.constants.f1;%fundamental resonance frequency
handles2.grho_val_min=(((handles2.drho_val_min./handles2.d2lam_val_max).*f1.*1.*cosd(handles2.phi_val_min./2)).^2)./1000;%calc. the grho value
set(handles2.grho_min,'string',handles2.grho_val_min,'ForegroundColor','r');
handles2.grho_min_txt.ForegroundColor='r';
handles2.d2lam_min.ForegroundColor='k';
handles2.d2lam_min_txt.ForegroundColor='k';
%maximum calculation
handles2.phi_val_max=str2double(handles2.phi_max.String);%phi (deg.)
handles2.drho_val_max=str2double(handles2.drho_max.String)./1000;%areal mass (kg/m^2)
handles2.d2lam_val_min=str2double(handles2.d2lam_min.String);
handles2.grho_val_max=(((handles2.drho_val_max./handles2.d2lam_val_min).*f1.*1.*cosd(handles2.phi_val_max./2)).^2)./1000;%calc. the grho value
set(handles2.grho_max,'string',handles2.grho_val_max,'ForegroundColor','r');
handles2.grho_max_txt.ForegroundColor='r';
handles2.d2lam_max.ForegroundColor='k';
handles2.d2lam_max_txt.ForegroundColor='k';

function d2lam_calc(~,~,handles,handles2)
%minimum calculation
handles2.phi_val_min=str2double(handles2.phi_min.String);%phi (deg.)
handles2.drho_val_min=str2double(handles2.drho_min.String)./1000;%areal mass (kg/m^2)
handles2.grho_val_max=str2double(handles2.grho_max.String).*1000;%grho values (Pa-kg/m^3)
f1=handles.din.constants.f1;%fundamental resonance frequency
handles2.d2lam_val_min=(handles2.drho_val_min.*f1.*1.*cosd(handles2.phi_val_min./2))./(sqrt(handles2.grho_val_max));%calc. the d2lam value
handles2.d2lam_min.String=handles2.d2lam_val_min;
handles2.d2lam_min.ForegroundColor='r';
handles2.d2lam_min_txt.ForegroundColor='r';
handles2.grho_min.ForegroundColor='k';
handles2.grho_min_txt.ForegroundColor='k';
%maximum calculation
handles2.phi_val_max=str2double(handles2.phi_max.String);%phi (deg.)
handles2.drho_val_max=str2double(handles2.drho_max.String)./1000;%areal mass (kg/m^2)
handles2.grho_val_min=str2double(handles2.grho_min.String).*1000;%grho values (Pa-kg/m^3)
handles2.d2lam_val_max=(handles2.drho_val_max.*f1.*1.*cosd(handles2.phi_val_max./2))./(sqrt(handles2.grho_val_min));%calc. the d2lam value
handles2.d2lam_max.String=handles2.d2lam_val_max;
handles2.d2lam_max.ForegroundColor='r';
handles2.d2lam_max_txt.ForegroundColor='r';
handles2.grho_max.ForegroundColor='k';
handles2.grho_max_txt.ForegroundColor='k';

function contour_Callback(hObject, ~, handles)
%This function creates contour plots of the master function.
set(handles.status,'string','Status: Generating contour plots...','backgroundcolor','k','foregroundcolor','r');
disp('Generating contour plots...');
harm_ratio_calc=handles.din.contour.harm_ratio_calc;
diss_ratio_calc=handles.din.contour.diss_ratio_calc;
f1=handles.din.constants.f1;%fundamental resonance freq
zq=handles.din.constants.zq;%quartz  load impedance

%viscoelastic functions
zqliq_n_calc=@(dfliq_n,dgliq_n) ((dfliq_n+1i.*dgliq_n).*pi.*zq)./(1i.*f1);%calculates the liquid load impedance based on the SLA eq.
drho_est=@(dfn,n) -zq.*dfn./(2.*n.*f1.^2);%estimated drho based on Sauerbrey equation
d2lam_n2_calc=@(d2lam_n1,n1,n2,phi) (d2lam_n1).*(n1./n2).^((phi./180)-1);%calculates the d2lambda ratio for harmonic n2
sauerbrey=@(n,drho) (2.*n.*(f1.^2).*drho)./(zq);%calculates the sauerbbrey frequency shift
zqliq_n2_calc=@(zqliq_n1,n1,n2) zqliq_n1.*sqrt(n2./n1);%calculates the liquid load impedance at harmonic n2
delfnstar_liq_calc=@(z_star_liq_n) (1i.*f1.*z_star_liq_n)./(pi.*zq);%calculates the complex frequency shift of a bare qcm xtal in liquid
Rliq_calc=@(delfnstar_liq,n,drho) (delfnstar_liq)./sauerbrey(n,drho);
Dn_calc=@(d2lam,phi) 2.*pi.*d2lam.*(1-1i.*tand(phi./2));
drho_calc=@(df,n,norm_delfstar) -drho_est(df,n)./real(norm_delfstar);%calculates the drho values
master2=@(Dn, Rliq) -((Dn.^-2)+(Rliq^2))./((cot(Dn)./Dn)+(Rliq));%master equation
rh_calc=@(d2lam_n1,phi,drho,n1,n2,delfnstar_liq1,delfnstar_liq2) real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n1,phi),phi),Rliq_calc(delfnstar_liq1,n1,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n2,phi),phi),Rliq_calc(delfnstar_liq2,n1,drho)));%harmonic ratio in terms of d2lam at n1
rd_calc=@(d2lam_n1,phi,drho,n1,n3,delfnstar_liq3) imag(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)));%dissipation ratio in terms of d2lam at n1
rh_calc2=@(d2lam_n1,phi,drho,n1,n2,delfnstar_liq1,delfnstar_liq2,nc) real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n1,phi),phi),Rliq_calc(delfnstar_liq1,n1,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n2,phi),phi),Rliq_calc(delfnstar_liq2,n1,drho)));%harmonic ratio in terms of designated d2lam
rd_calc2=@(d2lam_n1,phi,drho,n1,n3,delfnstar_liq3,nc) imag(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)));%dissipation ratio in terms of designated d2lam
grho_calc=@(drho,d2lam,f_n,phi) (drho.*f_n.*cosd(phi./2)./d2lam).^2;%in Pa-kg/m^3

%error functions
rh_var=@(covar1,dfn1,dfn2,n1,n2) ((n2./n1).^2./(dfn2.^2)).*(covar1(1,1).^2-2.*(dfn1./dfn2).*covar1(1,2)+(dfn1./dfn2).^2.*covar1(2,2));%variance of rh
rd_var=@(covar1,dfn1,dfn3,dgn3) ((dgn3./dfn3).^2.*covar1(3,3)-2.*(dgn3./dfn3).*covar1(3,4)+covar1(4,4))./dfn1.^2;%variance of rd
drho_var=@(covar1,n) (zq./(2.*n.*f1.^2)).^2.*covar1(1,1);%variance of drho
rh_rd_cov=@(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2)...
    (n2./n1)./(dfn2.*dfn3).*((covar1(1,3)-(dfn1./dfn2).*covar1(2,3)).*(-dgn3./dfn3)+covar1(1,4)-(dfn1./dfn2).*covar1(2,4));%covariance b/w rh and rd
rh_drho_cov=@(covar1,n,n1,n2,dfn1,dfn2) (zq./(2.*n.*f1.^2)).*(n2./n1).*(covar1(1,1)./dfn2-(dfn1./dfn2.^2).*covar1(1,2));%covariance between rh and drho
rd_drho_cov=@(covar1,dfn3,dgn3,n) (zq./(2.*n.*f1.^2))./dfn3.*(covar1(1,3).*(-dgn3./dfn3)+covar1(1,4));%covariance between rd and drho
ratio_cov=@(covar1,dfn1,dfn2,dfn3,dgn3,n,n1,n2) ...
    [rh_var(covar1,dfn1,dfn2,n1,n2), rh_rd_cov(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2), rh_drho_cov(covar1,n,n1,n2,dfn1,dfn2);...
    rh_rd_cov(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2), rd_var(covar1,dfn1,dfn3,dgn3), rd_drho_cov(covar1,dfn3,dgn3,n);...
    rh_drho_cov(covar1,n,n1,n2,dfn1,dfn2), rd_drho_cov(covar1,dfn3,dgn3,n), drho_var(covar1,n)];%covariance matrix of rh and rd (3x3 sym matrix)

%determine what n1, n2, n3 should be based on the value of the
%handles.contour_choice handle
[n1,n2,n3]=determine_harm_choice(get(handles.contour_choice,'value'));
%calculate experimental harmonic and dissipation ratio
df_n1=handles.din.cursor.(['interp_harmfi',num2str(n1)]);%Df_n1
df_n2=handles.din.cursor.(['interp_harmfi',num2str(n2)]);%Df_n2        
dg_n3=handles.din.cursor.(['interp_harmgi',num2str(n3)]);%Dg_n3
df_n3=handles.din.cursor.(['interp_harmfi',num2str(n3)]);%Df_n3
harm_ratio_exp=(n2./n1).*(df_n1./df_n2);%harmonic ratio from experiment
diss_ratio_exp=dg_n3./df_n3;%dissipation ratio from experiment
%calculate the error in the harm and diss ratios
try
    covariance_ratios=handles.din.stored_solutions.(['cov_ratios_',num2str(n1),num2str(n2),num2str(n3)]);%covariance matrix of the harm and diss ratio
catch
    covar1=handles.error_values.UserData;%Covariance matrix of raw freq and diss shifts
    covariance_ratios=ratio_cov(covar1,df_n1,df_n2,df_n3,dg_n3,n2,n1,n2);%covariance matrix of the harm and diss ratio
end
harm_ratio_error=2.*sqrt(covariance_ratios(1,1));
diss_ratio_error=2.*sqrt(covariance_ratios(2,2));
%calculate resonance frequencies at n1, n2, and n3
f_n1=f1*n1;%resonant frequency at n1 in Hz
f_n2=f1*n2;%resonant frequency at n2 in Hz
f_n3=f1*n3;%resonant frequency at n3 in Hz
guess_table=get(handles.uitable4,'data');
if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
    try
        drho=mean(unique([drho_est(df_n1,n1),drho_est(df_n2,n2),drho_est(df_n3,n3)]));%take the average value from all drho values and store as a guess value for drho (kg/m^2)
        if isnan(guess_values(3))==0             
            guess_table{1,2}=drho*1000;
            set(handles.uitable4,'data',guess_table);
        end%if isnan(drho_ave)==0
    catch
        drho=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
    end%try
else%if the manual guess for drho has been enabled
    drho=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
end%if get(handles.edit_drho,'value')==0
dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
n=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calculated
zqliq_n_calc=@(dfliq_n,dgliq_n) ((dfliq_n+1i.*dgliq_n).*pi.*zq)./(1i.*f1);%calculates the liquid load impedance based on the SLA eq.
zqliq_n2_calc=@(zqliq_n1,n1,n2) zqliq_n1.*sqrt(n2./n1);%calculates the liquid load impedance at harmonic n2
z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n);%calc the liq load impedance
z_star_liq_n1=zqliq_n2_calc(z_star_liq_n,n,n1);%liq load impedance at harmonic n1
z_star_liq_n2=zqliq_n2_calc(z_star_liq_n,n,n2);%liq load impedance at harmonic n2
z_star_liq_n3=zqliq_n2_calc(z_star_liq_n,n,n3);%liq load impedance at harmonic n3

if get(handles.d2lam_choice,'value')==1 %relative to d2lambda_n1
    nc=n1;
elseif get(handles.d2lam_choice,'value')==2 %relative to d2lambda_n2
    nc=n2;
end

if handles.qcm_map_sol.Value==1&&isfield(handles.din.stored_solutions,['n_',num2str(n1),'_',num2str(n1),num2str(n2),num2str(n3)])==1
    drho=handles.din.stored_solutions.(['n_',num2str(nc),'_',num2str(n1),num2str(n2),num2str(n3)]).drho(end,2)./1000;%calculate d2lamn value
    d2lamn=handles.din.stored_solutions.(['n_',num2str(nc),'_',num2str(n1),num2str(n2),num2str(n3)]).d2lamn(end,2);%calculate d2lamn value    
    d2lamn_error=handles.din.stored_solutions.(['n_',num2str(nc),'_',num2str(n1),num2str(n2),num2str(n3)]).error.d2lamn(end,2);%calculate d2lamn error
    phi=handles.din.stored_solutions.(['n_',num2str(nc),'_',num2str(n1),num2str(n2),num2str(n3)]).phi(end,2);%calculated phi value    
    phi_error=handles.din.stored_solutions.(['n_',num2str(nc),'_',num2str(n1),num2str(n2),num2str(n3)]).error.phi(end,2);%calculated phi error
    phi_var=linspace(phi-phi_error.*5,phi+phi_error.*5,200);
    d2lam_var=linspace(d2lamn-d2lamn_error.*5,d2lamn+d2lamn_error.*5,200);
    flag=1;
    cf_res1=linspace(harm_ratio_exp-harm_ratio_error*5,harm_ratio_exp+harm_ratio_error*5,200);
    cf_res2=linspace(diss_ratio_exp-diss_ratio_error*5,diss_ratio_exp+diss_ratio_error*5,200);
    
%     res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
%     phi_var=linspace(0,90,res);%define range of phi vales
%     d2lam_var=linspace(0,0.2,res);%define range of d2lam values
%     flag=0;
%     cf_res1=linspace(harm_ratio_exp-harm_ratio_error*500,harm_ratio_exp+harm_ratio_error*500,200);
%     cf_res2=linspace(diss_ratio_exp-diss_ratio_error*500,diss_ratio_exp+diss_ratio_error*500,200);
else
    res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
    phi_var=linspace(0,90,res);%define range of phi vales
    d2lam_var=linspace(0,0.2,res);%define range of d2lam values
    flag=0;
    cf_res1=linspace(harm_ratio_exp-harm_ratio_error*500,harm_ratio_exp+harm_ratio_error*500,200);
    cf_res2=linspace(diss_ratio_exp-diss_ratio_error*500,diss_ratio_exp+diss_ratio_error*500,200);
end

if handles.din.contour.label~=handles.din.guess_label%only regenerate the ratio matrices if the parameters have changed
    h=waitbar(0,'Generating plot...');%create progress bar
    harm_ratio_calc=nan(length(phi_var),length(d2lam_var));
    diss_ratio_calc=nan(length(phi_var),length(d2lam_var));
    for dum1=1:length(phi_var)
        try
            waitbar(dum1/length(phi_var),h);
        catch
            set(handles.status,'string','Cancelled plot generation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled plot generation.');
            return
        end%try
        for dum2=1:length(d2lam_var)       
            if get(handles.d2lam_choice,'value')==1 %relative to d2lambda_n1
                d2lam_n1=d2lam_var(dum2);%d2lambda at harmonic n1
                d2lam_n2=d2lam_n2_calc(d2lam_n1,n1,n2,phi_var(dum1));% calc. d2lambda ratio at harmonic n2
                d2lam_n3=d2lam_n2_calc(d2lam_n1,n1,n3,phi_var(dum1));% calc. d2lambda at harmonic n3      
                handles.din.n_ref=n1;
                harm_ratio_calc(dum1,dum2)=rh_calc2(d2lam_n1,phi_var(dum1),drho,n1,n2,delfnstar_liq_calc(z_star_liq_n1),delfnstar_liq_calc(z_star_liq_n2),n1);
                diss_ratio_calc(dum1,dum2)=rd_calc2(d2lam_n1,phi_var(dum1),drho,n1,n3,delfnstar_liq_calc(z_star_liq_n3),n1);
            elseif get(handles.d2lam_choice,'value')==2 %relative to d2lambda_n2              
                d2lam_n2=d2lam_var(dum2);% calc. d2lambada ratio at harmonic n2
                d2lam_n1=d2lam_n2_calc(d2lam_n2,n2,n1,phi_var(dum1));% d2lambda at harmonic n1
                d2lam_n3=d2lam_n2_calc(d2lam_n2,n2,n3,phi_var(dum1));% calc. d2lambda at harmonic n3
                handles.din.n_ref=n2;
                harm_ratio_calc(dum1,dum2)=rh_calc2(d2lam_n2,phi_var(dum1),drho,n1,n2,delfnstar_liq_calc(z_star_liq_n1),delfnstar_liq_calc(z_star_liq_n2),n2);
                diss_ratio_calc(dum1,dum2)=rd_calc2(d2lam_n2,phi_var(dum1),drho,n1,n3,delfnstar_liq_calc(z_star_liq_n3),n2);
            end% if get(handles.d2lam_choice,'value')==2
        end%for dum2=length(d2lam_var)
    end%for dum1=length(phi_var);
end%if handles.contour.label~=handles.din.guess_label

%store the calculated ratio in the handles structure so that the plot does
%not have to be reproduced every single time the callback function runs.
handles.din.contour.harm_ratio_calc=harm_ratio_calc;%store the harm_ratio_calc
handles.din.contour.diss_ratio_calc=diss_ratio_calc;%store the diss_ratio_calc
handles.din.contour.harm_ratio_exp=harm_ratio_exp;%store the harm_ratio_exp
handles.din.contour.harm_ratio_error=harm_ratio_error;%store the harm_ratio_error
handles.din.contour.diss_ratio_exp=diss_ratio_exp;%store the diss_ratio_exp
handles.din.contour.diss_ratio_error=diss_ratio_error;%store the diss_ratio_error
handles.din.contour.label=rand(1);%store a "label" for the stored calc. ratios
handles.din.guess_label=handles.din.contour.label;
guidata(handles.figure1,handles);%refresh the handles structure

%plot the contour map
clf(figure(1));
ff=figure(1);
set(ff,'inverthardcopy','off');
cmap=colormap;
if ff.Position(3)<=800
    set(ff,'position',[100+handles.figure1.Position(1) 400+handles.figure1.Position(2) 1120 420]);
end%if pos(3)<=0.71
s1=subplot(1,2,1);%plot the harmonic ratio map
contourf(d2lam_var,phi_var,harm_ratio_calc,cf_res1,'edgecolor','none');
xlabel(s1,['d/\lambda_',num2str(handles.din.n_ref)],'fontweight','bold','fontsize',14);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',14);
t1=title([num2str(n2),'\Deltaf_',num2str(n1),'/',num2str(n1),'\Deltaf_',num2str(n2)],'fontweight','bold','fontsize',14);
c1=colorbar; set(c1,'fontsize',12);
set(s1,'climmode','manual','color',cmap(1,:));

set(s1,'fontsize',14);
s2=subplot(1,2,2);%plot the dissipation ratio map
contourf(d2lam_var,phi_var,diss_ratio_calc,cf_res2,'edgecolor','none');
xlabel(s2,['d/\lambda_',num2str(handles.din.n_ref)],'fontweight','bold','fontsize',14);
ylabel(s2,'\phi (deg.)','fontweight','bold','fontsize',14);
t2=title(['\Delta\Gamma_',num2str(n3),'/\Deltaf_',num2str(n3)],'fontweight','bold','fontsize',14);
set(s2,'fontsize',14);
linkaxes([s1,s2],'xy');%link the axes
set(handles.status,'string','Status: Plot succesfully generated!','foregroundcolor','r','backgroundcolor','k');
c2=colorbar; set(c2,'fontsize',12);
set(s2,'climmode','manual','color',cmap(1,:));

%plot the solutions
if get(handles.qcm_map_sol,'value')==1
    hold(s1,'on');
    [map1,hc1]=contour(s1,d2lam_var,phi_var,harm_ratio_calc,[harm_ratio_exp-harm_ratio_error harm_ratio_exp+harm_ratio_error],'edgecolor',[.3 .3 .3],'linestyle','-');
    [map2,hc2]=contour(s1,d2lam_var,phi_var,diss_ratio_calc,[diss_ratio_exp-diss_ratio_error diss_ratio_exp+diss_ratio_error],'edgecolor',[.3 .3 .3],'linestyle','--');
    clabel(map1,hc1,'fontweight','bold');
    clabel(map2,hc2,'fontweight','bold');    
    hold(s2,'on');
    [map3,hc3]=contour(s2,d2lam_var,phi_var,harm_ratio_calc,[harm_ratio_exp-harm_ratio_error harm_ratio_exp+harm_ratio_error],'edgecolor',[.3 .3 .3],'linestyle','-');
    [map4,hc4]=contour(s2,d2lam_var,phi_var,diss_ratio_calc,[diss_ratio_exp-diss_ratio_error diss_ratio_exp+diss_ratio_error],'edgecolor',[.3 .3 .3],'linestyle','--');
    clabel(map3,hc3,'fontweight','bold');
    clabel(map4,hc4,'fontweight','bold');    
    t1.String=[t1.String,' = ',num2str(harm_ratio_exp),'\pm',num2str(harm_ratio_error)];
    t2.String=[t2.String,' = ',num2str(diss_ratio_exp),'\pm',num2str(diss_ratio_error)];
    clabel(map1,hc1,'labelspacing',200);
    clabel(map2,hc2,'labelspacing',200);
    clabel(map3,hc3,'labelspacing',200);
    clabel(map4,hc4,'labelspacing',200);
    if flag==1
        errorbar(s1,d2lamn,phi,phi_error,phi_error,d2lamn_error,d2lamn_error,'mo','linewidth',1.5,'markersize',8);
        errorbar(s2,d2lamn,phi,phi_error,phi_error,d2lamn_error,d2lamn_error,'mo','linewidth',1.5,'markersize',8);
    end%flag==1
end%if get(handles.qcm_map_sol,'value')==1
linkaxes([s1 s2],'xy');
drawnow
try
    delete(h);
end
disp('Plot succesfully generated');



function d2lam_choice_Callback(hObject, ~, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
guidata(hObject,handles);

function contour_choice_Callback(hObject, ~, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
guidata(hObject,handles);

function plot_solutions_Callback(hObject, ~, handles)
%this function plots out the calculted solutions as a function of time
%NOTE this program currently only plots out data calculated at n1
set(handles.status,'string','Status: Plotting solutions...','backgroundcolor','k','foregroundcolor','r');
disp('Plotting solutions...');

%plot viscoelastic parameters
ff1=figure(4);clf(ff1);
set(ff1,'position',[handles.figure1.Position(1),handles.figure1.Position(2),1466,414]);
s1=subplot(1,3,1);%plot drho
s2=subplot(1,3,2);%plot grho
s3=subplot(1,3,3);%plot phi
[L_string,xstr]=plot_viscoelastic(s1,s2,s3,handles,handles.din.n1);%plot the viscoelastic properties at the first harmonic (5MHz)
xlabel(s1,xstr,'fontweight','bold','fontsize',14);
ylabel(s1,'d\rho (g/m^2)','fontweight','bold','fontsize',14);
xlabel(s2,xstr,'fontweight','bold','fontsize',14);
ylabel(s2,['|G^*_',num2str(handles.din.n1),'|\rho (Pa-g/cm^3)'],'fontweight','bold','fontsize',14);
xlabel(s3,xstr,'fontweight','bold','fontsize',14);
ylabel(s3,['\phi_',num2str(handles.din.n1),' (deg.)'],'fontweight','bold','fontsize',14);
set(s1,'fontsize',14);
set(s2,'fontsize',14);
set(s3,'fontsize',14);
L=legend(s3,'','','','','','');
set(L,'string',L_string);

%plot the raw frequencies
ff2=figure(3);clf(ff2);
set(ff2,'position',[handles.figure1.Position(1),500+handles.figure1.Position(2),1466,414]);
r1=subplot(1,2,1);
r2=subplot(1,2,2);
[L0_string,~]=plot_raw_freq(r1,r2,handles);
%plot the calc. frequencies
[L_string,xstr]=plot_calc_freq(r1,r2,handles,0);
[~,~]=plot_calc_freq(r1,r2,handles,1);
%format the plot
L=legend(r2,'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','');
set(L,'string',[L0_string,L_string]);
xlabel(r1,xstr,'fontweight','bold','fontsize',14);
ylabel(r1,'\Deltaf (Hz)','fontweight','bold','fontsize',14);
xlabel(r2,xstr,'fontweight','bold','fontsize',14);
ylabel(r2,'\Delta\Gamma (Hz)','fontweight','bold','fontsize',14);
set(r1,'fontsize',14);
set(r2,'fontsize',14);
if get(handles.dg_log,'value')==1;set(r2,'yscale','log');else set(r2,'yscale','linear');end;
if get(handles.linkx,'value')==1%if axes linking is turned on
    linkaxes([r1,r2,s1,s2,s3],'x');
end%if get(handles.linkx,'value')==1
%store the error values in Df and Dg in the figures
data=get(handles.error_values,'userdata');
data1=cell(6,3);
for dum=1:6
    data1(dum,:)=[{dum*2-1},{data(dum,1)},{data(dum,2)}];
end
data1=[{'harmonic order'},{'df error (Hz)'},{'dG error (Hz)'};data1];
set(figure(3),'userdata',data1);
set(figure(4),'userdata',data1);
set(handles.status,'string','Status: Solutions plotted!','backgroundcolor','k','foregroundcolor','r');
disp('Solutions plotted!');


function error_var=thin(d2lam,phi)
Dn=@(d2lam,phi) 2.*pi.*d2lam.*(1-1i.*tand(phi./2));
master=@(d2lam,phi) -tan(Dn(d2lam,phi))./Dn(d2lam,phi);
thin=@(d2lam,phi) -(1+(1./3).*Dn(d2lam,phi).^2);
error=@(d2lam,phi) ((thin(d2lam,phi)-master(d2lam,phi))./(master(d2lam,phi))).*100;
rd1=@(d2lam,phi) imag(master(d2lam,phi))./real(master(d2lam,phi));
rd2=@(d2lam,phi) imag(thin(d2lam,phi))./real(thin(d2lam,phi));
error2=@(d2lam,phi) ((rd2(d2lam,phi)-rd1(d2lam,phi))./rd1(d2lam,phi)).*100;
Denolf=@(d2lam,phi) -(0.229).*(d2lam.^2).*phi;%from 2011 paper
error3=@(d2lam,phi) ((Denolf(d2lam,phi)-rd1(d2lam,phi))./rd1(d2lam,phi)).*100;

res=length(d2lam);
error_var=nan(res,res);
error_var2=nan(res,res);
Denolf_var=nan(res,res);
thin_var=nan(res,res);
for dum=1:length(d2lam)
    for dum2=1:length(phi)
        thin_var(dum2,dum)=thin(d2lam(dum),phi(dum2));
        error_var(dum2,dum)=error(d2lam(dum),phi(dum2));
        error_var2(dum2,dum)=error2(d2lam(dum),phi(dum2));
        Denolf_var(dum2,dum)=error3(d2lam(dum),phi(dum2));
    end
end

function qcm_map_Callback(hObject, eventdata, handles)
%this function plots out the contour plots for the master equation
set(handles.status,'string','Status: Generating QCM maps...','backgroundcolor','k','foregroundcolor','r');
disp('Generating QCM maps...');
qcm_map=handles.din.qcm_map;
res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
phi_var=linspace(0,90,res);%define range of phi vales
d2lam_var=linspace(0,1,res);%define range of d2lam values
f1=handles.din.constants.f1;%fundamental resonance freq
zq=handles.din.constants.zq;%quartz  load impedance
f_n=f1*str2double(get(handles.dgliq_harm,'string'));
guess_table=get(handles.uitable4,'data');

if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
    try
        drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
        drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
        drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
        drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values        
        guess_table{1,2}=drho_ave;
        set(handles.din.uitable4,'data',guess_table);
    catch
        drho_ave=guess_table{1,2}./1000;%kg/m^2
    end
else%if the manual guess for drho has been enabled
drho_ave=guess_table{1,2}./1000;%kg/m^2
end%if get(handles.edit_drho,'value')==0
dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
n=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at  
z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance
h=waitbar(0,'Generating plot...');%create progress bar    
if handles.din.qcm_label~=handles.din.guess_label2%recalculate the QCM map matrix if the parameters were changed    
    %calculate the normalized complex freq. shifts
    for dum1=1:length(phi_var);
        try waitbar(dum1/(length(phi_var)+10),h);
        catch
            set(handles.status,'string','Cancelled plot generation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled plot generation.');
            return
        end%try
        for dum2=1:length(d2lam_var)
            qcm_map(dum1,dum2)=master(d2lam_var(dum2),phi_var(dum1),z_star_liq_n,f_n,drho_ave);            
        end%for dum2=1:length(d2lam_var)
    end%for dum1=1:length(phi_var);
end%if handles.din.qcm_label~=handles.din.guess_label2

%store and update variables in the handles structure
handles.din.qcm_map=qcm_map;
handles.din.qcm_label=rand(1);
handles.din.guess_label2=handles.din.qcm_label;
guidata(hObject,handles);
%plot the QCM map
clf(figure(2));
ff=figure(2);
pos=get(ff,'position');
if pos(3)<=800
    set(ff,'position',[pos(1)*.5 pos(2) pos(3)*2 pos(4)]);
end%if pos(3)<=0.71
s1=subplot(1,2,1);%plot the real component
contourf(d2lam_var,phi_var,real(qcm_map),[-logspace(-2,2,500),logspace(-2,0,500)],'edgecolor','none');hold on;
decay2lambda=@(x) cotd(x./2)./(2.*pi);
plot(s1,decay2lambda(phi_var),phi_var,'k--');
[~,temp]=contour(s1,d2lam_var,phi_var,real(thin(d2lam_var,phi_var)),[-5 5],'edgecolor','k','linestyle','--','showtext','on');
set(s1,'xlim',[0 d2lam_var(end)]);
xlabel(s1,'d/\lambda_n','fontweight','bold','fontsize',12);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Deltaf_n','/\Deltaf_s_n'],'fontweight','bold','fontsize',12);
colormap(parula); colorbar;
waitbar(.95,h);
s2=subplot(1,2,2);%plot the imaginary component
contourf(d2lam_var,phi_var,imag(qcm_map),logspace(-10,0,500),'edgecolor','none');hold on;
plot(s2,decay2lambda(phi_var),phi_var,'k--');
[~,temp]=contour(s2,d2lam_var,phi_var,imag(thin(d2lam_var,phi_var)),[0 5],'edgecolor','k','linestyle','--','showtext','on');
set(s2,'xlim',[0 d2lam_var(end)]);
xlabel(s2,'d/\lambda_n','fontweight','bold','fontsize',12);
ylabel(s2,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Delta\Gamma_n','/\Deltaf_s_n'],'fontweight','bold','fontsize',12);
cmap=colormap(parula);

%%THIS CODE PLOTS OUT RELATIVE ERROR OF THIN FILM APPROXIMATION (3RD ORDER
%%TAYLOR EXPANSION ABOUT D2LAM=0) FOR +/-5% RELATIVE TO THE MASTER EQUATION
% f3.f=figure(3);clf(f3.f);
% f3.s1=subplot(1,2,1);hold(f3.s1,'on');
% f3.s2=subplot(1,2,2);hold(f3.s2,'on');
% contourf(f3.s1,d2lam_var,phi_var,real(thin(d2lam_var,phi_var)),linspace(-100,100,200),'edgecolor','none');
% v=caxis(f3.s1);
% [C,temp]=contour(f3.s1,d2lam_var,phi_var,real(thin(d2lam_var,phi_var)),'edgecolor','r','linestyle','--');
% temp.LevelList=[-5 5];
% [~,temp]=contour(f3.s1,d2lam_var,phi_var,imag(thin(d2lam_var,phi_var)),'edgecolor','k','linestyle','--');
% temp.LevelList=[-5 5];
% caxis(f3.s1,v);
% contourf(f3.s2,d2lam_var,phi_var,imag(thin(d2lam_var,phi_var)),linspace(-100,100,200),'edgecolor','none');
% v=caxis(f3.s2);
% [~,temp]=contour(f3.s2,d2lam_var,phi_var,real(thin(d2lam_var,phi_var)),'edgecolor','r','linestyle','--');
% temp.LevelList=[-5 5];
% [~,temp]=contour(f3.s2,d2lam_var,phi_var,imag(thin(d2lam_var,phi_var)),'edgecolor','k','linestyle','--');
% temp.LevelList=[-5 5];
% caxis(f3.s2,v);
% xlabel(f3.s1,'d/\lambda_n','fontweight','bold','fontsize',12);
% ylabel(f3.s1,'\phi (deg.)','fontweight','bold','fontsize',12);
% xlabel(f3.s2,'d/\lambda_n','fontweight','bold','fontsize',12);
% ylabel(f3.s2,'\phi (deg.)','fontweight','bold','fontsize',12);
% colorbar(f3.s1);
% colorbar(f3.s2);

% cmap=[ones(1000,3).*cmap(1,:);cmap];
% colormap(s1,cmap);
caxis(s1,[-2 0.5])
caxis(s2,[0 1]);
colorbar;
waitbar(1,h);
cf=z_star_liq_n/(drho_ave*f_n);
set(ff,'name',['drho=',num2str(drho_ave.*1000),'g/m^2  cf=',num2str(cf)]);
set(handles.status,'string','Status: Plot succesfully generated!','foregroundcolor','r','backgroundcolor','k');
if get(handles.qcm_map_sol,'value')==1%plot solutions if the option to do so is turned on
    [L_string]=plot_d2lam_phi(handles,handles.din.n1,s1,s2);
    L=legend(s2,'','','','','','','','','','','','','location','best');
    set(L,'string',[{''},L_string])
end%if get(handles.qcm_map_sol,'value')==1
try  delete(h); end%try
disp('Plot succesfully generated');

function [L_string]=plot_d2lam_phi(handles,n1,s1,s2)
%this function plots the viscoelastic parrameters in subplots s1, s2, and
%s3 at harmonic n1
sols=fieldnames(handles.din.stored_solutions);%extract out the filenames (n_<harmonic at which the parameters were calcaulated at>_<harmonics used to calc data>
count=1;
L_string=[];
for dum=1:length(sols)
    test=sols{dum};
    test1=strfind(test,'_');
    index=str2double(sols{dum}(test1(1)+1:test1(end)-1));
    if index==n1
        index2(count)=dum;
        count=count+1;
    end%if index==1
end%for dum=1:length(sols)
for dum=1:length(index2)
    temp=sols{index2(dum)};
    ydata1=handles.din.stored_solutions.(temp).d2lamn(:,2);%d2lam at n=1
    ydata2=handles.din.stored_solutions.(temp).phi(:,2);%phi at n=1
    edata1=handles.din.stored_solutions.(temp).error.d2lam(:,2);% the error for drho
    edata2=handles.din.stored_solutions.(temp).error.phi(:,2);% the error for phi
    k=strfind(temp,'_');
    L_string=[L_string,{temp(k(end)+1:end)}];%string for legend box
    try
        errorbar(s1,ydata1,ydata2,edata2,edata2,edata1,edata1,...
            'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        errorbar(s2,ydata1,ydata2,edata2,edata2,edata1,edata1,...
            'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
    catch
        set(handles.status,'string',['Status: Could not plot errorbar dataset for ',temp(end-2:end)],'backgroundcolor','y','foregroundcolor','r');
        disp(['Could not plot errorbar dataset for ',temp(end-2:end)]);
        plot(s1,ydata1,ydata2,...
            'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        plot(s2,ydata1,ydata2,...
            'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
    end%try        
end%for dum=1:length(index2)

function [L_string,xstr]=plot_viscoelastic(s1,s2,s3,handles,n1)
%this function plots the viscoelastic parrameters in subplots s1, s2, and
%s3 at harmonic n1
sols=fieldnames(handles.din.stored_solutions);%extract out the filenames (n_<harmonic at which the parameters were calcaulated at>_<harmonics used to calc data>
count=1;
L_string=[];
for dum=1:length(sols)
    test=sols{dum};
    test1=strfind(test,'_');
    index=str2double(sols{dum}(test1(1)+1:test1(end)-1));
    if index==n1
        index2(count)=dum;
        count=count+1;
    end%if index==1
end%for dum=1:length(sols)
for dum=1:length(index2)
    temp=sols{index2(dum)};
    [xstr,xdata]=determine_time_units(handles.din.stored_solutions.(temp).drho(:,1),handles);
    ydata1=handles.din.stored_solutions.(temp).drho(:,2);%drho at n=1
    ydata2=handles.din.stored_solutions.(temp).grho(:,2);%grho at n=1
    ydata3=handles.din.stored_solutions.(temp).phi(:,2);%phi at n=1
    edata1=handles.din.stored_solutions.(temp).error.drho(:,2);
    edata2=handles.din.stored_solutions.(temp).error.grho(:,2);
    edata3=handles.din.stored_solutions.(temp).error.phi(:,2);
    r1=find(ydata1~=0); ydata1=ydata1(r1); edata1=edata1(r1);
    r2=find(ydata2~=0); ydata2=ydata2(r2); edata2=edata2(r2);
    r3=find(ydata3~=0); ydata3=ydata3(r3); edata3=edata3(r3);
    xdata=xdata(r1);
    if get(handles.norm_time,'value')==1%run this code if the option to normalize the timepoints is turned on
        xdata=xdata-min(xdata);
    end%if get(handles.norm_time,'value')==1
    k=strfind(temp,'_');
    L_string=[L_string,{temp(k(end)+1:end)}];%string for legend box
    try
        errorbar(s1,xdata,ydata1,2.*edata1,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        errorbar(s2,xdata,ydata2,2.*edata2,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
        errorbar(s3,xdata,ydata3,2.*edata3,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s3,'on');
    catch
        set(handles.status,'string',['Status: Could not plot errorbar dataset for ',temp(end-2:end)],'backgroundcolor','y','foregroundcolor','r');
        disp(['Could not plot errorbar dataset for ',temp(end-2:end)]);
        plot(s1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        plot(s2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
        plot(s3,xdata,ydata3,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s3,'on');
    end%try        
end%for dum=1:length(index2)
set(s1,'ylim',[min(ydata1)-abs(min(ydata1)*.15) max(ydata1)+abs(min(ydata1)*.15)]);
set(s2,'ylim',[min(ydata2)-abs(min(ydata2)*.5) max(ydata2)+abs(min(ydata2)*.5)]);
set(s3,'ylim',[min(ydata3)-abs(min(ydata3)*.5) max(ydata3)+abs(min(ydata3)*.5)]);


function [L0_string,xstr]=plot_raw_freq(r1,r2,handles)
%this function plots the raw frequences in subplots r1 and r2
L0_string=[];
freq_shifts=handles.din.rawdata.freq_shift;
for dum=1:length(handles.din.harmtot)
    [xstr,xdata]=determine_time_units(freq_shifts(:,1),handles);%determine xdata
    if get(handles.norm_time,'value')==1%run this code if the option to normalize the timepoints is turned on
        xdata=xdata-xdata(1);
    end%  if get(handles.norm_time,'value')==1
    plot(r1,xdata(:,1),freq_shifts(:,handles.din.harmtot(dum)+1),'-','color',handles.din.marker_color{0.5*(handles.din.harmtot(dum)+1)},'linewidth',1);hold(r1,'on');%plot Df
    plot(r2,xdata(:,1),freq_shifts(:,handles.din.harmtot(dum)+2),'-','color',handles.din.marker_color{0.5*(handles.din.harmtot(dum)+1)},'linewidth',1);hold(r2,'on');%plot pred. Dg
    L0_string=[L0_string,{['n=',num2str(handles.din.harmtot(dum))]}];
end%for dum=1:length(handles.din.harmtot)


function [L_string,xstr]=plot_calc_freq(r1,r2,handles,flag)
%this function plots the calculated frequencies in subplots r1 and r2
L_string=[];
marker_style=handles.din.marker_style;
marker_color=handles.din.marker_color;
sols=fieldnames(handles.din.stored_solutions);%extract out the filenames (n_<harmonic at which the parameters were calcaulated at>_<harmonics used to calc data>
for dum1=1:2:11%calculate for each harmonic
    count=1;
    index2=[];
    for dum=1:length(sols)
        test=sols{dum};
        test1=strfind(test,'_');
        index=str2double(sols{dum}(test1(1)+1:test1(end)-1));
            if index==dum1
                index2(count)=dum;
                count=count+1;
            end%if index==1
    end%for dum=1:length(sols)
    for dum=1:length(index2)%calc for each dataset combination
        temp=sols{index2(dum)};
        test1=strfind(temp,'_');
        [xstr,xdata]=determine_time_units(handles.din.stored_solutions.(temp).g_shifts(:,1),handles);%determine xdata
        if get(handles.norm_time,'value')==1%run this cord if the option to normalize the timepoints is turned on
            xdata=xdata-min(xdata);
        end%  if get(handles.norm_time,'value')==1
        ydata1=handles.din.stored_solutions.(temp).f_shifts(:,3);%calc Df
%         ldata1=ydata1-real(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,2));
%         udata1=real(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,3))-ydata1;
        ydata2=handles.din.stored_solutions.(temp).g_shifts(:,3);%Calc Dg       
%         ldata2=ydata2-imag(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,2));
%         udata2=imag(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,3))-ydata2;
        k=strfind(temp,'_');
        %get rid of 0,0 datapoints
        non1=find(xdata==0);
        if isempty(non1)==0
            for bum=1:length(non1)
                if ydata1(non1(bum))==0||ydata2(non1(bum))==0
                    xdata(non1(bum))=nan;
                    ydata1(non1(bum))=nan;
                    ydata2(non1(bum))=nan;
                    ldata1(non1(bum))=nan;
                    ldata2(non1(bum))=nan;
                    udata1(non1(bum))=nan;
                    udata2(non1(bum))=nan;
                end%if ydata1(non1(bum))==0                
            end%for bum=1:length(non1)
        end%if isempty(non1)==0
        if dum1==1
            L_string=[L_string,{[temp(k(end)+1:end)]}];%string for legend box
        end%if dum1==1
        if flag==0&&dum1==1%if the flag is turned off, plot the pred freq black
            hold(r1,'on');hold(r1,'on');
            plot(r1,xdata(1),ydata1(1),'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none');
            plot(r2,xdata(1),ydata2(1),'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none'); 
        elseif flag==1
            index3=0.5*(dum1+1);
            hold(r1,'on');hold(r1,'on');
%             errorbar(r1,xdata,ydata1,ldata1,udata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
%             errorbar(r2,xdata,ydata2,ldata2,udata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none'); 
            plot(r1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
            plot(r2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none'); 
        end%if flag==0
    end%for dum=1:length(index2)
end%for dum1=1:2:11

function qcm_map_sol_Callback(hObject, eventdata, handles)
handles.din.contour.label=rand(1);%store a "label" for the stored calc. ratios
guidata(handles.figure1,handles);%refresh the handles structure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACK FUNCTIONS IN THE ACTIVE HARMONICS PANEL

function harm1_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);
function harm3_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);
function harm5_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);
function harm7_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);
function harm9_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);
function harm11_Callback(hObject, ~, handles)
myharmfcn(hObject,handles);

function myharmfcn(hObject,handles)
set(handles.status,'string','Loading harmonic...');
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
set(handles.status,'string','Status: Ready!');

function harmtot=active_harm(handles)
%this function determines which harmonics are being actively used, which is
%determined based on the toggle state of the radio handles. in the
%active harmonics panel
harmtot=[];
for dum=1:2:11
    name=['harm',num2str(dum)];
    if get(handles.(name),'value')==1
        harmtot=[harmtot;dum];
    end%if get(handles.(name),'value');==1
end%for dum=1:2:11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACK FUNCTIONS IN THE PLOTTING PANEL
function handles=plot_fcn(handles,rawdata)
%this function is a customized plotting function that plots the rawdata
%that was imported by the user.
cla(handles.axes1);cla(handles.axes2);
legend_str=[];
for dum=1:length(handles.din.harmtot)
    namef=['harmf',num2str(handles.din.harmtot(dum))];
    nameg=['harmg',num2str(handles.din.harmtot(dum))];
    handles.din.(namef)=rawdata.freq_shift(:,handles.din.harmtot(dum)+1);
    handles.din.(nameg)=rawdata.freq_shift(:,handles.din.harmtot(dum)+2);
    I_f=find(isnan(handles.din.(namef))==0);%exclude elements that are NaN
    I_g=find(isnan(handles.din.(nameg))==0);%exclude elements that are NaN        
    [xstring,time]=determine_time_units(rawdata.freq_shift(:,1),handles);%calculate time array
    handles.din.(namef)=[time(I_f),handles.din.(namef)(I_f)];%time (col 1) and delta f (col 2) for harmonic dum with NaN excluded
    handles.din.(nameg)=[time(I_g),handles.din.(nameg)(I_g)];%time (col 1) and delta gamma (col 2) for harmonic dum with NaN excluded
    if get(handles.norm_time,'value')==0%determine whether or not to normalize te time scale
        hold (handles.axes1,'on');%hold handles.axes1 on the "on" state
        plot(handles.axes1,handles.din.(namef)(:,1),handles.din.(namef)(:,2)./handles.din.harmtot(dum),...
            'o-','markersize',10,'linewidth',2,'color',handles.din.marker_color{(handles.din.harmtot(dum)+1)/2});%Delta_f versus time
        xlabel(handles.axes1,xstring,'fontweight','bold','fontsize',12);
        ylabel(handles.axes1,'\Deltaf/n (Hz)','fontweight','bold','fontsize',12);
        hold(handles.axes2,'on');%hold handles.axes2 on the "off" state
        plot(handles.axes2,handles.din.(nameg)(:,1),handles.din.(nameg)(:,2),...
            'o-','markersize',10,'linewidth',2,'color',handles.din.marker_color{(handles.din.harmtot(dum)+1)/2});%Delta_gamma versus time
        xlabel(handles.axes2,xstring,'fontweight','bold','fontsize',12);
        ylabel(handles.axes2,'\Delta\Gamma (Hz)','fontweight','bold','fontsize',12);
        legend_str=[legend_str,{['n=',num2str(handles.din.harmtot(dum))]}];%write a string array that will be used for the legend box
    else%plot with normalized time scale (initial timepoint will be set to 0)
        hold (handles.axes1,'on');%hold handles.axes1 on the "on" state
        plot(handles.axes1,handles.din.(namef)(:,1)-min(handles.din.(namef)(:,1)),handles.din.(namef)(:,2)./handles.din.harmtot(dum),...
            'o-','markersize',10,'linewidth',2,'color',handles.din.marker_color{(handles.din.harmtot(dum)+1)/2});%Delta_f versus time
        xlabel(handles.axes1,xstring,'fontweight','bold','fontsize',12);
        ylabel(handles.axes1,'\Deltaf/n (Hz)','fontweight','bold','fontsize',12);
        hold(handles.axes2,'on');%hold handles.axes2 on the "off" state
        plot(handles.axes2,handles.din.(nameg)(:,1)-min(handles.din.(namef)(:,1)),handles.din.(nameg)(:,2),...
            'o-','markersize',10,'linewidth',2,'color',handles.din.marker_color{(handles.din.harmtot(dum)+1)/2});%Delta_gamma versus time
        xlabel(handles.axes2,xstring,'fontweight','bold','fontsize',12);
        ylabel(handles.axes2,'\Delta\Gamma (Hz)','fontweight','bold','fontsize',12);
        legend_str=[legend_str,{['n=',num2str(handles.din.harmtot(dum))]}];%write a string array that will be used for the legend box
    end%if get(handles.norm_time,'value')==0    
end%for dum=1:length(handles.din.harmtot)
%draw legend entries
L=legend(handles.axes1,'','','','','','');%allocate 6 empty entries (note that the warning message that normally pops up has been suppressed)
set(L,'units','normalized','string',legend_str,'location','southeast');
L=legend(handles.axes2,'','','','','','');%allocate 6 empty entries (note that the warning message that normally pops up has been suppressed)
set(L,'units','normalized','string',legend_str,'location','southeast');
%determine whether or not to link the x axis
if get(handles.linkx,'value')==1
    linkaxes([handles.axes1,handles.axes2],'x');
else
    linkaxes([handles.axes1,handles.axes2],'off');
end%if get(handles.linkx,'value')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACKS ASSOCIATED WITH THE BARE CRYSTAL READING
%PANEL
function handles=bare_crystal(handles)
rawdata=handles.din.rawdata;%extract the imported rawdata from handles structure
freq_shift_ref=rawdata.freq_shift_ref';%extract out the reference frequency shifts
set(handles.bare_table,'data',freq_shift_ref);%export reference data into the table in the bare crystal reading panel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS ALL OF THE CALLBACK FUNCTIONS ASSOCIATED WITH THE
%BUTTONS ON THE GUI TOOLBAR
function debug_ClickedCallback(hObject, eventdata, handles)
%This function allows the user to access the handles structure. To eit out
%of the debugging mode, click on "Exit debuging mode" in the top ribbon or
%type in "return" in the command window. See the MATLAB keyboard function
%for more details.
keyboard;


function home_button_ClickedCallback(hObject, eventdata, handles)
%This function sets the gui figure back to its original state.
handles=home_state(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function bare_xtal_ClickedCallback(hObject, eventdata, handles)
%This function allows the user to load bare_crystal data and recalculate
%the frequency shifts relative to the bare_crystal data.
disp('Importing bare .mat file');
set(handles.status,'string','Status: Importing bare crystal .mat file');
[filename,pathfile,~]=uigetfile('.mat', 'Load bare crystal datafile',handles.din.bare_path);%prompt user to choose the bare crystal file
set(handles.status,'string','Status: Importing .mat file...');drawnow;
if isempty(filename)%run this code if the user cancels out of loading the bare crystal file
    set(handles.status,'string','Status: Unable to load bare crystal datafile',...
        'foregroundcolor','k','backgroundcolor','r');
    return
end%if isempty(filename)
try
    load([pathfile,filename]);%load the datafile
    for dum1=2:size(abs_freq,2)%find the number of columns in the abs_freq variable loaded from the bare crystal file and calculate the average
        temp=abs_freq(:,dum1);%extract out the designated column
        if mod(dum1,2)==0&&isnan(nanmean(temp))==0%if the column is associated with delta f
            handles.din.rawdata.ref_freq(dum1/2)=nanmean(temp);%calculate the average and redefine the reference harmonic frequencies
        elseif mod(dum1,2)==1&&isnan(nanmean(temp))==0%if the column is associated with delta gamma
            handles.din.rawdata.ref_diss((dum1-1)/2)=nanmean(temp);%calculate the average and redefine the reference harmonic dissipations
        end%if mod(dum1,2)==0&&isnan(mean(temp))==0%if the column is associated with delta f
    end%for dum1=2:size(abs_freq,2)
    %update the bare_crystal reading table
    set(handles.bare_table,'data',[handles.din.rawdata.ref_freq',handles.din.rawdata.ref_diss']);
    disp('Reference frequency and dissipation values have been redefined.');
    disp(['ref_freq: ',num2str(handles.din.rawdata.ref_freq)]);
    disp(['ref_diss: ',num2str(handles.din.rawdata.ref_diss)]);
    %redefine the frequency shifts relative to the loaded bare crystal reference values
    handles.din.rawdata.freq_shift=[handles.din.rawdata.freq_shift(:,1),...
        handles.din.rawdata.abs_freq(:,2)-handles.din.rawdata.ref_freq(1),handles.din.rawdata.abs_freq(:,3)-handles.din.rawdata.ref_diss(:,1),...
        handles.din.rawdata.abs_freq(:,4)-handles.din.rawdata.ref_freq(2),handles.din.rawdata.abs_freq(:,5)-handles.din.rawdata.ref_diss(:,2),...
        handles.din.rawdata.abs_freq(:,6)-handles.din.rawdata.ref_freq(3),handles.din.rawdata.abs_freq(:,7)-handles.din.rawdata.ref_diss(:,3),...
        handles.din.rawdata.abs_freq(:,8)-handles.din.rawdata.ref_freq(4),handles.din.rawdata.abs_freq(:,9)-handles.din.rawdata.ref_diss(:,4),...
        handles.din.rawdata.abs_freq(:,10)-handles.din.rawdata.ref_freq(5),handles.din.rawdata.abs_freq(:,11)-handles.din.rawdata.ref_diss(:,5),...
        handles.din.rawdata.abs_freq(:,12)-handles.din.rawdata.ref_freq(6),handles.din.rawdata.abs_freq(:,13)-handles.din.rawdata.ref_diss(:,6)];
    handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
    set(handles.status,'string',['Status: Bare crystal file, ',filename,', successfully loaded!']);
    disp(['Status: Bare crystal file, ',filename,', successfully loaded!']);
    handles.din.bare_flag=1;%turn on the flag, showing that a bare crystal was loaded
    handles.din.bare_path=pathfile;
    set(hObject,'tooltipstring',['<html>Bare crystal loaded on ',datestr(clock),'<br />Filename: ',filename,'<br />Filepath: ',pathfile,'</html>']);
    
    guidata(hObject, handles);
catch err_message
    assignin('base','err_message',err_message);
    set(handles.status,'string','Status: Error in loading bare crystal file!');
    disp('Error in loading bare crystal file!');
    return
end%try



% --------------------------------------------------------------------
function immersed_bare_xtal_ClickedCallback(hObject, eventdata, handles)
% This fcn accepts an immersed bare xtal file (bare xtal immersed in a
% semi-infinite layer). An average is taken for each harmonic and the
% values are set for the dissipation shift associated with the approprate
% harmonic for the semi-infinite medium.
disp('Importing immersed bare .mat file');
set(handles.status,'string','Status: Importing immersed bare crystal .mat file');
[filename,pathfile,~]=uigetfile('.mat', 'Load bare crystal datafile',handles.din.bare_path);%prompt user to choose the bare crystal file
set(handles.status,'string','Status: Importing .mat file...');drawnow;
if isempty(filename)%run this code if the user cancels out of loading the bare crystal file
    set(handles.status,'string','Status: Unable to load bare crystal datafile',...
        'foregroundcolor','k','backgroundcolor','r');
    return
end%if isempty(filename)
try
    load([pathfile,filename]);%load the datafile
    for dum1=2:size(abs_freq,2)%find the number of columns in the abs_freq variable loaded from the bare crystal file and calculate the average
        temp=abs_freq(:,dum1);%extract out the designated column
        if mod(dum1,2)==0&&isnan(nanmean(temp))==0%if the column is associated with delta f
            handles.din.immersed.immersed_freq(dum1/2)=nanmean(temp);%calculate the average and store the immersed harmonic frequencies
        elseif mod(dum1,2)==1&&isnan(nanmean(temp))==0%if the column is associated with delta gamma
            handles.din.immersed.immersed_diss((dum1-1)/2)=nanmean(temp);%calculate the average and store th immersed harmonic dissipations
        end%if mod(dum1,2)==0&&isnan(mean(temp))==0%if the column is associated with delta f
    end%for dum1=2:size(abs_freq,2)
    index=(str2double(get(handles.dgliq_harm,'string'))+1)/2;
    set(handles.dgliq,'string',handles.din.immersed.immersed_diss(index));
catch err_message
    assignin('base','err_message',err_message);
    set(handles.status,'string','Status: Error in loading immersed bare crystal file!');
    disp('Error in loading immersed bare crystal file!');
    return
end%try
guidata(handles.figure1,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION DEALS WITH CALLBACK FUNCTIONS ASSOCIATED WITH THE CALCULATION
%OPTIONS PANEL

function radio1_Callback(hObject, ~, handles)
function radio2_Callback(hObject, ~, handles)
function radio3_Callback(hObject, ~, handles)
function radio4_Callback(hObject, ~, handles)
function radio5_Callback(hObject, ~, handles)
function radio6_Callback(hObject, ~, handles)
function harmchoice1_Callback(hObject, ~, handles)
function harmchoice2_Callback(hObject, ~, handles)
function harmchoice3_Callback(hObject, ~, handles)
function harmchoice4_Callback(hObject, ~, handles)
function harmchoice5_Callback(hObject, ~, handles)
function harmchoice6_Callback(hObject, ~, handles)

function initial_time_Callback(hObject, ~, handles)
timepoints=handles.din.rawdata.freq_shift(:,1);%xtract out the timepoints
ind=find(isnan(timepoints)==0);
timepoints=timepoints(ind);
handles.din.extent=length(timepoints);
ind=find(isnan(timepoints)==0);
timepoints=timepoints(ind);
[~,timepoints]=determine_time_units(timepoints,handles);
ind1=str2double(get(hObject,'string'));
if get(handles.norm_time,'value')==0
    set(hObject,'tooltipstring',['<html>Initial timepoint: ',num2str(timepoints(ind1)),' min.','<br/>'...
        'Total # of pts: ',num2str(length(timepoints)),'</html>']);
else
    set(hObject,'tooltipstring',['<html>Initial timepoint: ',num2str(timepoints(ind1)-min(timepoints)),' min.','<br/>',...
        'Total # of pts: ',num2str(length(timepoints)),'</html>']);
end%if get(handles.norm_time,'value')==0

function extent_Callback(hObject, ~, handles)
timepoints=handles.din.rawdata.freq_shift(:,1);%xtract out the timepoints
ind=find(isnan(timepoints)==0);
timepoints=timepoints(ind);
handles.din.extent=length(timepoints);
ind=find(isnan(timepoints)==0);
timepoints=timepoints(ind);
[~,timepoints]=determine_time_units(timepoints,handles);
ind1=str2double(get(hObject,'string'));
if get(handles.norm_time,'value')==0
    set(hObject,'tooltipstring',['<html>Final timepoint: ',num2str(timepoints(ind1)),' min.','<br/>'...
        'Total # of pts: ',num2str(length(timepoints)),'</html>']);
else
    set(hObject,'tooltipstring',['<html>Final timepoint: ',num2str(timepoints(ind1)-min(timepoints)),' min.','<br/>',...
        'Total # of pts: ',num2str(length(timepoints)),'</html>']);
end%if get(handles.norm_time,'value')==0

function dgliq_Callback(hObject, ~, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters
handles.din.guess_label2=rand(1);%"relabel" the guess parameters
guidata(hObject,handles);

function dgliq_harm_Callback(hObject, ~, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
handles.din.guess_label2=rand(1);%"relabel" the guess parameters
index=(str2double(get(hObject,'string'))+1)/2;
set(handles.dgliq,'string',handles.din.immersed.immersed_diss(index));
guidata(hObject,handles);

function edit_drho_Callback(hObject, ~, handles)
guess_table=get(handles.uitable4,'data');
if get(hObject,'value')==1%this signals to the program to not auto calculate the drho guess, but to use the user input value    
    guess_table{1,1}='<HTML><body text="0000FF"><b>d&rho; (g/m<sup>2</sup>)</b></body></html>';    
elseif get(hObject,'value')==0
    guess_table{1,1}='<HTML><body text="000000"><b>d&rho; (g/m<sup>2</sup>)</b></body></html>';    
end%if get(hObject,'value')==1
set(handles.uitable4,'data',guess_table);

function cursor_Callback(hObject, ~, handles)
%This function allows the user to pick a point from the plots. The
%datapoint is then used for analysis.
handles.din.stored_solutions=[];%clear out the stored solutions
button=1;%this is a flag that represent the which button was pushed during datapoint selection
count=1;
handles.din.harmtot=active_harm(handles);
guidata(handles.figure1,handles);%immediately update the handles structure
set(handles.status,'string','Status: Click R mouse button o close out of cursor mode. Click L mouse button to continue selection','backgroundcolor','k','foregroundcolor','r');
disp('Click R mouse button to close out of cursor mode. Click L mouse button to continue selection');
while button==1
    try
        figure(handles.figure1);%make the GUI figure the current figure
        [selected_time0,~,button]=ginput(1);%this turns the pointer to a cursor, allowing the user to select a series of datapoints from the axes plots
        if isempty(selected_time0)==1
            set(handles.status,'string','Status: Selected timepoint is empty, please select another timepoint','backgroundcolor','y','foregroundcolor','r');
            disp('Selected timepoint is empty, please select another timepoint');
        end%if isempty(selected_time0)==1
        switch get(handles.time_units,'value')%convert the timepoint unit to minutes
            case 2 %if the xaxis is in units of hours
                selected_time0=selected_time0.*60;
            case 3% if the xaxis is in units of days
                selected_time0=selected_time0.*(60.*24);
        end%switch xaxis_units
        if get(handles.norm_time,'value')==1%run this code if the time scale has been normalized
            selected_time0=min(handles.din.rawdata.freq_shift(:,1))+selected_time0;%convert the relative timepoint to an absolute time point
        end%if get(handles.norm_time,'value')==1
        for dum=1:length(handles.din.harmtot)% run this code for each active harmonics that is shown in the axes plots
            namefi=['interp_harmfi',num2str(handles.din.harmtot(dum))];
            namegi=['interp_harmgi',num2str(handles.din.harmtot(dum))];
            xdata_f=handles.din.rawdata.freq_shift(:,1);%rawdata time array associated with harmonic harmtot(dum)
            ydata_f=handles.din.rawdata.freq_shift(:,handles.din.harmtot(dum)+1);%rawdata \Delta f array associated with harmonic harmtot(dum)
            xdata_g=handles.din.rawdata.freq_shift(:,1);%rawdata time array associated with harmonic harmtot(dum)
            ydata_g=handles.din.rawdata.freq_shift(:,handles.din.harmtot(dum)+2);%rawdata \Delta g array associated with harmonic harmtot(dum)
            %get rid of nans
            ind=find(isnan(xdata_f)==0);
            xdata_f=xdata_f(ind);
            ydata_f=ydata_f(ind);
            xdata_g=xdata_g(ind);
            ydata_g=ydata_g(ind);
            ind2=find(isnan(ydata_f)==0);
            xdata_f=xdata_f(ind2);
            ydata_f=ydata_f(ind2);
            xdata_g=xdata_g(ind2);
            ydata_g=ydata_g(ind2);            
            xaxis_units=get(handles.time_units,'value');%determine what units the timepoint is in                        
            try%sometimes interp1 does not always work (especially if the dataset is not monotonic)
                handles.din.cursor.(namefi)=interp1(xdata_f,ydata_f,selected_time0,'linear');%interpolate the frequency shift defined by selected_time0
                handles.din.cursor.(namegi)=interp1(xdata_g,ydata_g,selected_time0,'linear');%interpolate the frequency shift defined by selected_time0
            catch%if it does not work try a different method (find nearest timepoint);
                ind3=find(abs(xdata_f-selected_time0)==min(abs(xdata_f-selected_time0)),1,'first');%find the index associated with the nearest timepoint
                handles.din.cursor.(namefi)=ydata_f(ind3);%find the frequency of the nearest point relative to the query point
                handles.din.cursor.(namegi)=ydata_g(ind3);%find the dissipation fo the nearest point relative to the query point
                selected_time0=xdata_g(ind3);%redefine the query timepoint to the timepoint that is nearest to the query timepoint
            end%try
            handles.din.cursor.selected_time=selected_time0;%store the selected timepoint in handles structure
            pred_freq((handles.din.harmtot(dum)+1)/2,:)=[{handles.din.cursor.(namefi)},{['']},{handles.din.cursor.(namegi)},{['']}];%output the interpolated frequency shifts
            set(handles.uitable3,'data',pred_freq);%display the interpolated freq shifts in the table
            set(handles.uitable2,'data',[]);%clear the table in the viscoelastic parameters panel
            %export the selected and interpolated datapoints onto the base
            %workspace (for debugging purposes)
            assignin('base',namefi,handles.din.cursor.(namefi));
            assignin('base',namegi,handles.din.cursor.(namegi));
            assignin('base','selected_time0',selected_time0);
        end%for dum=1:length(handles.din.harmtot)
        if get(handles.auto_solve,'value')==1
            handles=qcm_solver(handles,count);%solve for viscoelastic parameters            
        end%if get(handles.auto_solve,'value')==1
        %--------------------insert contour solutions script here
        figure(handles.figure1);%make the GUI figure the current figure
        count=count+1;
    catch err_msg
        assignin('base','err_msg',err_msg);
        disp('Error with data selection!');
        set(handles.status,'string','Status: Error with data selection!','backgroundcolor','r','foregroundcolor','k');
    end%try
end%while button==1
handles.din.contour.label=rand(1);
guidata(hObject,handles);
set(handles.status,'string','Status: Exited out of cursor mode. Ready...','backgroundcolor','k','foregroundcolor','r');
disp('Exited out of cursor mode. Ready...');

function [radiotot,harmchoice]=radio_total(handles)
%this function calculates the number and type of dataset to be used for
%calculation
radiotot=[];
harmchoice=[];
for dum=1:6
    name=['radio',num2str(dum)];
    name2=['harmchoice',num2str(dum)];
    if get(handles.(name),'value')==1
        radiotot=[radiotot;dum];
        harmchoice=[harmchoice;get(handles.(name2),'value')];
    end%if get(handles.(name),'value')==1
end%for dum=1:6

function solve_all_Callback(hObject, ~, handles)
%this function solves every other datapoint based on the increment defined
%by solve_inc
handles.din.stored_solutions=[];%clear out the stored solutions
set(handles.status,'string','Status: Solving...','backgroundcolor','k','foregroundcolor','r');
disp('Solving...');
raw_freq_shifts=handles.din.rawdata.freq_shift;%extract out the freq shift data
timepoints=raw_freq_shifts(:,1);%xtract out the timepoints
%figure out the total number of datapoints
It=find(isnan(timepoints)==0);
timepoints=timepoints(It);
raw_freq_shifts=raw_freq_shifts(It,:);
count=1;
h=waitbar(0,'Calculating');%create status bar
initial=str2double(get(handles.initial_time,'string'));%determine the start time index
extent=str2double(get(handles.extent,'string'));%determine the end time index
for dum=initial:str2double(get(handles.solve_inc,'string')):extent%run this code for each timepoint
    try
        waitbar(dum/extent,h,['Calculating...',num2str(dum),' of ',num2str(extent)]);
    catch
        set(handles.status,'string','Status: Calculation cancelled!','backgroundcolor','k','foregroundcolor','r');
        disp('Calculation cancelled!');
        return
    end%try
    for dum2=1:length(handles.din.harmtot)
        try
            namefi=['interp_harmfi',num2str(handles.din.harmtot(dum2))];
            namegi=['interp_harmgi',num2str(handles.din.harmtot(dum2))];
            handles.din.cursor.selected_time=timepoints(dum);%extractout the timepoint
            handles.din.cursor.(namefi)=raw_freq_shifts(dum,handles.din.harmtot(dum2)+1);%extract the associated Df
            handles.din.cursor.(namegi)=raw_freq_shifts(dum,handles.din.harmtot(dum2)+2);%extract the associated Dg
            pred_freq((handles.din.harmtot(dum2)+1)/2,:)=[{handles.din.cursor.(namefi)},{['']},{handles.din.cursor.(namegi)},{['']}];%output the interpolated frequency shifts
            set(handles.uitable3,'data',pred_freq);%display the interpolated freq shifts in the table
            set(handles.uitable2,'data',[]);%clear the table in the viscoelastic parameters panel
            if isnan(handles.din.cursor.(namefi))==1||isnan(handles.din.cursor.(namegi))==1
                disp(['Frequency shifts at timepoint: ',num2str(handles.din.cursor.selected_time),' min. are NaN']);
                break
            else
            end%if isnan(handles.din.cursor.(namefi))==1||isnan(handles.din.cursor.(namegi))==1      
        catch err_msg
            assignin('base','err_msg',err_msg);
            set(handles.status,'string','Status: Unable to find solution!','backgroundcolor','k','foregroundcolor','r');
            disp('Unable to find solution!');
        end%try
        guidata(hObject,handles);
    end%dum2=1:length(handles.din.harmtot)    
    handles=qcm_solver(handles,count);%solve for viscoelastic parameters  
    count=count+1;
end%for dum=length(timepoints)
try delete(h); end%try to delete the status bar
handles.din.contour.label=rand(1);%store a "label" for the stored calc. ratios
guidata(hObject,handles);
set(handles.status,'string','Status: Calculation finished!','backgroundcolor','k','foregroundcolor','r');
disp('Calculation finished!');

function solve_inc_Callback(hObject, eventdata, handles)


function auto_solve_Callback(hObject, eventdata, handles)


function handles=grho_guess_Callback(hObject, ~, handles)
%Running this callback function will calculate the predicted d2lambda value.
data=get(handles.uitable4,'data');%extract out values from guess table
phi=cell2mat(data(2,2));%phi (deg.)
drho=cell2mat(data(1,2))./1000;%areal mass, drho (kg/m^2)
grho=cell2mat(data(3,2)).*1000;%grho values (Pa-kg/m^3)
f1=handles.din.constants.f1;%fundamental resonance frequency
d2lam=(drho.*f1.*1.*cosd(phi./2))./(sqrt(grho));%calc. the d2lam value
data{4,1}='<HTML><body text=#FF0000><b>d/&lambda;</b></body></html>';
data{3,1}='<HTML><body text=#000000><b>&#8739;G*&#8739;&rho; (Pa-g/cm<sup><font size=2>3</font></sup>)</b></body></html>';
data{4,2}=d2lam;
set(handles.uitable4,'data',data);
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
handles.din.guess_label2=rand(1);
guidata(hObject,handles);

function handles=d2lam_guess_Callback(hObject, ~, handles)
%Running this callback function will calculate the predicted Grho value.
data=get(handles.uitable4,'data');%extract out values from guess table
d2lam=cell2mat(data(4,2));%d/lambda value
phi=cell2mat(data(2,2));%phi (deg.)
drho=cell2mat(data(1,2))./1000;%areal mass, drho (kg/m^2)
f1=handles.din.constants.f1;%fundamental resonance frequency
grho=(((drho./d2lam).*f1.*1.*cosd(phi./2)).^2)./1000;%calc. the grho value
data{3,1}='<HTML><body text=#FF0000><b>&#8739;G*&#8739;&rho; (Pa-g/cm<sup><font size=2>3</font></sup>)</b></body></html>';
data{4,1}='<HTML><body text=#000000><b>d/&lambda;</b></body></html>';
data{3,2}=grho;
set(handles.uitable4,'data',data);
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
handles.din.guess_label2=rand(1);
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QCM EQUATIONS

function handles=qcm_solver(handles,count)
%this function calculates the two variables, d2lambda and phi.
%handles: the gui handles structure
%count: keeps track which row to store the calculated viscoelasetic
%parameters in a variable
[radiotot,harmchoice]=radio_total(handles);%determine how many datasets and which harmonics that will be used to perform the calculation
f1=handles.din.constants.f1;%fundamentalr resonance freq in Hz
zq=handles.din.constants.zq;% the quartz load impedance

%viscoelastic functions
zqliq_n_calc=@(dfliq_n,dgliq_n) ((dfliq_n+1i.*dgliq_n).*pi.*zq)./(1i.*f1);%calculates the liquid load impedance based on the SLA eq.
drho_est=@(dfn,n) -zq.*dfn./(2.*n.*f1.^2);%estimated drho based on Sauerbrey equation
d2lam_n2_calc=@(d2lam_n1,n1,n2,phi) (d2lam_n1).*(n1./n2).^((phi./180)-1);%calculates the d2lambda ratio for harmonic n2
sauerbrey=@(n,drho) (2.*n.*(f1.^2).*drho)./(zq);%calculates the sauerbbrey frequency shift
zqliq_n2_calc=@(zqliq_n1,n1,n2) zqliq_n1.*sqrt(n2./n1);%calculates the liquid load impedance at harmonic n2
delfnstar_liq_calc=@(z_star_liq_n) (1i.*f1.*z_star_liq_n)./(pi.*zq);%calculates the complex frequency shift of a bare qcm xtal in liquid
Rliq_calc=@(delfnstar_liq,n,drho) (delfnstar_liq)./sauerbrey(n,drho);
Dn_calc=@(d2lam,phi) 2.*pi.*d2lam.*(1-1i.*tand(phi./2));
drho_calc=@(df,n,norm_delfstar) -drho_est(df,n)./real(norm_delfstar);%calculates the drho values
master2=@(Dn, Rliq) -((Dn.^-2)+(Rliq^2))./((cot(Dn)./Dn)+(Rliq));%master equation
rh_calc=@(d2lam_n1,phi,drho,n1,n2,delfnstar_liq1,delfnstar_liq2) real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n1,phi),phi),Rliq_calc(delfnstar_liq1,n1,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n2,phi),phi),Rliq_calc(delfnstar_liq2,n1,drho)));%harmonic ratio in terms of d2lam at n1
rd_calc=@(d2lam_n1,phi,drho,n1,n3,delfnstar_liq3) imag(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,n1,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)));%dissipation ratio in terms of d2lam at n1
rh_calc2=@(d2lam_n1,phi,drho,n1,n2,delfnstar_liq1,delfnstar_liq2,nc) real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n1,phi),phi),Rliq_calc(delfnstar_liq1,n1,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n2,phi),phi),Rliq_calc(delfnstar_liq2,n1,drho)));%harmonic ratio in terms of designated d2lam
rd_calc2=@(d2lam_n1,phi,drho,n1,n3,delfnstar_liq3,nc) imag(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)))./...
    real(master2(Dn_calc(d2lam_n2_calc(d2lam_n1,nc,n3,phi),phi),Rliq_calc(delfnstar_liq3,n3,drho)));%dissipation ratio in terms of designated d2lam
grho_calc=@(drho,d2lam,f_n,phi) (drho.*f_n.*cosd(phi./2)./d2lam).^2;%in Pa-kg/m^3

%error functions
rh_var=@(covar1,dfn1,dfn2,n1,n2) ((n2./n1).^2./(dfn2.^2)).*(covar1(1,1).^2-2.*(dfn1./dfn2).*covar1(1,2)+(dfn1./dfn2).^2.*covar1(2,2));%variance of rh
rd_var=@(covar1,dfn1,dfn3,dgn3) ((dgn3./dfn3).^2.*covar1(3,3)-2.*(dgn3./dfn3).*covar1(3,4)+covar1(4,4))./dfn1.^2;%variance of rd
drho_var=@(covar1,n) (zq./(2.*n.*f1.^2)).^2.*covar1(1,1);%variance of drho
rh_rd_cov=@(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2)...
    (n2./n1)./(dfn2.*dfn3).*((covar1(1,3)-(dfn1./dfn2).*covar1(2,3)).*(-dgn3./dfn3)+covar1(1,4)-(dfn1./dfn2).*covar1(2,4));%covariance b/w rh and rd
rh_drho_cov=@(covar1,n,n1,n2,dfn1,dfn2) (zq./(2.*n.*f1.^2)).*(n2./n1).*(covar1(1,1)./dfn2-(dfn1./dfn2.^2).*covar1(1,2));%covariance between rh and drho
rd_drho_cov=@(covar1,dfn3,dgn3,n) (zq./(2.*n.*f1.^2))./dfn3.*(covar1(1,3).*(-dgn3./dfn3)+covar1(1,4));%covariance between rd and drho
ratio_cov=@(covar1,dfn1,dfn2,dfn3,dgn3,n,n1,n2) ...
    [rh_var(covar1,dfn1,dfn2,n1,n2), rh_rd_cov(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2), rh_drho_cov(covar1,n,n1,n2,dfn1,dfn2);...
    rh_rd_cov(covar1,dfn1,dfn2,dfn3,dgn3,n1,n2), rd_var(covar1,dfn1,dfn3,dgn3), rd_drho_cov(covar1,dfn3,dgn3,n);...
    rh_drho_cov(covar1,n,n1,n2,dfn1,dfn2), rd_drho_cov(covar1,dfn3,dgn3,n), drho_var(covar1,n)];%covariance matrix of rh and rd (3x3 sym matrix)

%define anonymous functions to determine the error in d2lam_n
J11_d2lam_n=@(n,n1,phi) (n/n1).^(1-phi./180);%derivative of d2lam_n with respect to d2lam_n1 (where n1 is the harmonic in which the solution was determined at)
J12_d2lam_n=@(n,n1,phi,d2lamn1) 2./pi.*(-1/180).*d2lamn1.*log(n/n1).*(n/n1).^(1-phi.*180);%derivative of d2lam_n with respect to phi
var_d2lam=@(n,n1,phi,d2lamn1,covar_x) J11_d2lam_n(n,n1,phi).^2.*covar_x(1,1)+...
    2.*J11_d2lam_n(n,n1,phi).*J12_d2lam_n(n,n1,phi,d2lamn1).*covar_x(1,2)+...
    J12_d2lam_n(n,n1,phi,d2lamn1).^2.*covar_x(2,2);%variance of d2lam_n

%define anonymous dunctions to determine the error in grho
J11_grho_n1=@(n1,phi,d2lam_n1,drho) -2.*(drho.*n1.*f1.*cosd(phi./2)).^2.*d2lam_n1.^(-3);%derivative of grho with respect to d2lam_n1
J12_grho_n1=@(n1,phi,d2lam_n1,drho) -0.5.*drho.^2.*(n1.*f1).^2.*sind(phi).*d2lam_n1.^(-2).*pi./180;%derivative of grho with respect to phi
J13_grho_n1=@(n1,phi,d2lam_n1,drho) 2.*drho.*(n1.*f1).^2.*cosd(phi./2).^2.*d2lam_n1.^(-2);%derivative of grho with respect to drho
var_grho_n1=@(n1,phi,d2lam_n1,drho,covar_x) J11_grho_n1(n1,phi,d2lam_n1,drho).^2.*covar_x(1,1)+...
    J12_grho_n1(n1,phi,d2lam_n1,drho).^2.*covar_x(2,2)+J13_grho_n1(n1,phi,d2lam_n1,drho).^2.*covar_x(3,3)+...
    2.*J11_grho_n1(n1,phi,d2lam_n1,drho).*J12_grho_n1(n1,phi,d2lam_n1,drho).*covar_x(1,2)+...
    2.*J12_grho_n1(n1,phi,d2lam_n1,drho).*J13_grho_n1(n1,phi,d2lam_n1,drho).*covar_x(2,3)+...
    2.*J11_grho_n1(n1,phi,d2lam_n1,drho).*J13_grho_n1(n1,phi,d2lam_n1,drho).*covar_x(1,3);%variance of grho at n1 (the harmonic in which the solution was determined at)

for dum0=1:length(radiotot)
    calc_table=nan(6,9);    %preallocate the calculation table
    handles.din.solved.harmchoice=harmchoice(dum0);
    [n1,n2,n3]=determine_harm_choice(harmchoice(dum0));%determine what harmonic datasets will be used to calculate viscoelastic parameters
    df_n1=handles.din.cursor.(['interp_harmfi',num2str(n1)]);%Df_n1
    df_n2=handles.din.cursor.(['interp_harmfi',num2str(n2)]);%Df_n2        
    dg_n3=handles.din.cursor.(['interp_harmgi',num2str(n3)]);%Dg_n3
    df_n3=handles.din.cursor.(['interp_harmfi',num2str(n3)]);%Df_n3
    
    %calculate the liquid load impedance
    dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
    dfliq_n=-dgliq_n;%freq shift used for Newtonian liq load impedance calc
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at            
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n);%calc the liq load impedance at dgliq_harm
    
    %extract out the guess parameters
    guess_table=get(handles.uitable4,'data');
    d2lam_n1=guess_table{4,2};%guess value for d2lambda at n1
    phi=guess_table{2,2};%guess value for phi
    if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
        try
            drho=mean(unique([drho_est(df_n1,n1),drho_est(df_n2,n2),drho_est(df_n3,n3)]));%take the average value from all drho values and store as a guess value for drho (kg/m^2)
            if isnan(guess_values(3))==0             
                guess_table{1,2}=drho*1000;
                set(handles.uitable4,'data',guess_table);
            end%if isnan(drho_ave)==0
        catch
            drho=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
        end%try
    else%if the manual guess for drho has been enabled
        drho=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
    end%if get(handles.edit_drho,'value')==0
    if isnan(drho)==1%check to see if drho is nan, if so, revert back to default value
        drho=1;
    end%if if isnan(drho)==1
    
    %numerically solve for the variables
    try        
        temp_cov=handles.error_values.UserData;%Covariance matrix of raw freq and diss shifts
        covar1=[temp_cov(n1,n1) temp_cov(n1,n2) temp_cov(n1,n3) temp_cov(n1,n3+1);...
            temp_cov(n2,n1) temp_cov(n2,n2) temp_cov(n2,n3) temp_cov(n2,n3+1);...
            temp_cov(n3,n1) temp_cov(n3,n2) temp_cov(n3,n3) temp_cov(n3,n3+1);...
            temp_cov(n3+1,n1) temp_cov(n3+1,n2) temp_cov(n3+1,n3) temp_cov(n3+1,n3+1)];
        harm_ratio_exp=(n2./n1).*(df_n1./df_n2);%harmonic ratio from experiment
        diss_ratio_exp=dg_n3./df_n3;%dissipation ratio from experiment                        
        covariance_ratios=ratio_cov(covar1,df_n1,df_n2,df_n3,dg_n3,n2,n1,n2);%covariance matrix of the harm and diss ratio
        handles.din.stored_solutions.(['cov_ratios_',num2str(n1),num2str(n2),num2str(n3)])=covariance_ratios;%harmonic and dissipation ratio covariance matrix
        
        %solve for phi and d2lambda
        guess_values=[d2lam_n1,phi,drho];            
        fcns=@(x) [(rh_calc2(x(1),x(2),x(3),n1,n2,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)),delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n2)),n1)-harm_ratio_exp);...
            (rd_calc2(x(1),x(2),x(3),n1,n3,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n3)),n1)-diss_ratio_exp);...
            (1000.*(drho_calc(df_n1,n1,master2(Dn_calc(x(1),x(2)),Rliq_calc(delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)),n1,x(3))))-x(3)))];%calculates the difference between the calculated and experimental ratios
        
        handles.din.fit_options=optimset('Display','off','tolfun',1e-8,'tolx',1e-8,'Maxiter',1e8,'MaxFunEvals',1e8);
        [solved,resnorm,residual,exitflag,output,~,jacobian]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,handles.din.fit_options);%numerically solved
        output.resnorm=resnorm;% stores the squared 2-norm of the residue at x
        output.residual=residual;%stores the value of the residual at the solution x
        output.jacobian=jacobian;%store the jacobian matrix from fsolve
        output.exitflag=exitflag;%store the exitflag that was outputed from fsolve
        disp(['exitflag: ',num2str(exitflag),', resnorm: ',num2str(resnorm),', # of pts solved: ',num2str(count),', dataset: ',num2str(n1),num2str(n2),num2str(n3)]);
        disp(solved);
        if exitflag~=1&&resnorm>0.0001%stop the code if the exitflag is not 1 or if the sum of the squares of the residuals are greater than 1e-10
            set(handles.status,'string','Status: Exitflag is not 1!','backgroundcolor','y','foregroundcolor','r');
            disp('Exitflag is not 1!');
            continue
        else
            set(handles.status,'string','Calculating...','backgroundcolor','k','foregroundcolor','r');
        end%if exitflag~=1
        handles.din.solved.output_misc_from_fsolve=output;%store the output structure in the handles structure
        handles.din.solved.solution=solved;%store the solution, solved(1): d2lam at n1 and solved(2): phi, solved(3): drho at n1        
            
        %solve for the other viscoelastic parameters                            
            unique_n=unique([n1;n2;n3;active_harm(handles)],'stable');
        for dum=1:length(unique_n)
            name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
            handles.din.solved.(name).d2lam=d2lam_n2_calc(solved(1),n1,unique_n(dum),solved(2));%store the d2lambda ratio at n1
            handles.din.solved.(name).norm_delfstar=master2(Dn_calc(handles.din.solved.(name).d2lam,solved(2)),Rliq_calc(delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,unique_n(dum))),unique_n(dum),solved(3)));
            handles.din.solved.(name).drho=solved(3).*1000;
            handles.din.solved.(name).complex_pred_freq=handles.din.solved.(name).norm_delfstar.*sauerbrey(unique_n(dum),solved(3));%calculate the complex frequency shifts (Hz)
            handles.din.solved.(name).phi=solved(2);%store the viscoelastic phase angle at n1 (degrees)
            handles.din.solved.(name).grho=grho_calc(handles.din.solved.(name).drho./1000,handles.din.solved.(name).d2lam,unique_n(dum)*f1,solved(2))./1000;%calculate grho at n1 in units of Pa-g/cm^3
            handles.din.solved.(name).cf=zqliq_n2_calc(z_star_liq_n,dgliq_harm,unique_n(dum))./(unique_n(dum).*f1.*handles.din.solved.(name).drho);%calculate the "cf" term
            handles.din.solved.(name).lam_rho=(sqrt(handles.din.solved.(name).grho.*1000)./(unique_n(dum).*f1.*cosd(solved(2)/2))).*1000;%calculate \lambda\rho (um-g/cm^3)
            handles.din.solved.(name).del2d=-1/imag(handles.din.solved.(name).d2lam.*2.*pi.*(1-1i.*tand(handles.din.solved.(name).phi./2)));
            handles.din.solved.(name).rh_calc=rh_calc2(solved(1),solved(2),solved(3),n1,n2,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)),delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n2)),n1);
            handles.din.solved.(name).rd_calc=rd_calc2(solved(1),solved(2),solved(3),n1,n3,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n3)),n1);
            handles.din.solved.(name).z_star_liq_n=z_star_liq_n;%liquid acoustic impedance
            handles.din.solved.(name).dgliq_harm=dgliq_harm;
            handles.din.solved.(name).delfnstar_liq_n1=delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1));
            handles.din.solved.(name).delfnstar_liq_n2=delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n2));
            handles.din.solved.(name).delfnstar_liq_n3=delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n3));                       
            
            %store the calculated values in a matrix that will be used to output the results in the handles.uitable2
            calc_table((unique_n(dum)+1)/2,:)=[handles.din.solved.(name).drho,handles.din.solved.(name).grho,handles.din.solved.(name).phi,...
                handles.din.solved.(name).d2lam,handles.din.solved.(name).lam_rho,real(handles.din.solved.(name).norm_delfstar),...
                imag(handles.din.solved.(name).norm_delfstar),abs(handles.din.solved.(name).cf),handles.din.solved.(name).del2d];
            %store the predicted freq shifts in a matrix that will be used to output the results in the handles.uitable3
            pred_shifts((unique_n(dum)+1)/2,:)=[handles.din.cursor.(['interp_harmfi',num2str(unique_n(dum))]),real(handles.din.solved.(name).complex_pred_freq),...
                handles.din.cursor.(['interp_harmgi',num2str(unique_n(dum))]),imag(handles.din.solved.(name).complex_pred_freq)];
            
            handles.din.stored_solutions.(name).d2lamn(count,:)=[handles.din.cursor.selected_time,handles.din.solved.(name).d2lam];
            handles.din.stored_solutions.(name).drho(count,:)=[handles.din.cursor.selected_time,handles.din.solved.(name).drho];%store the drho value and timepoint
            handles.din.stored_solutions.(name).grho(count,:)=[handles.din.cursor.selected_time,handles.din.solved.(name).grho];%store the grho value and timepoint
            handles.din.stored_solutions.(name).phi(count,:)=[handles.din.cursor.selected_time,handles.din.solved.(name).phi];%store the phi value and timepoint
            handles.din.stored_solutions.(name).f_shifts(count,:)=[handles.din.cursor.selected_time,pred_shifts((unique_n(dum)+1)/2,1:2)];%store the experimental (col2) and predicted freq shifts (col3) and timepoint (col1: timepoint, col2: exp_shift, col3: calc_shift)
            handles.din.stored_solutions.(name).g_shifts(count,:)=[handles.din.cursor.selected_time,pred_shifts((unique_n(dum)+1)/2,3:4)];%store the experimental (col2) and predicted freq shifts (col3) and timepoint (col1: timepoint, col2: exp_shifts, col3: clac shifts)                        
        end%for dum=1:length(unique_n)
        %update the data tables in the gui
        set(handles.uitable2,'data',calc_table);%update the viscoelastic parameter table
        set(handles.uitable3,'data',pred_shifts);%update the predicted frequency shift table                        
    catch err_msg
        assignin('base','err_msg',err_msg);
        set(handles.status,'string','Status: Error in solving!','backgroundcolor','r','foregroundcolor','k');
        disp('Error in solving!');
        return
    end%try
    
    try%calculate the error
        covar_x1=full(output.jacobian)\covariance_ratios/transpose(full(output.jacobian));%determine the covariance of numerically solved solution vector @ d2lam @ n1
        disp(' ');
        for dum=1:length(unique_n)
            name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
            fcns=@(x) [rh_calc2(x(1),x(2),x(3),n1,n2,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)),delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n2)),unique_n(dum))-harm_ratio_exp;...
                rd_calc2(x(1),x(2),x(3),n1,n3,delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n3)),unique_n(dum))-diss_ratio_exp;...
                1000.*(drho_calc(df_n1,n1,master2(Dn_calc(x(1),x(2)),Rliq_calc(delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)),n1,x(3))))-x(3))];%calculates the difference between the calculated and experimental ratios                        
                
            [solved2,~,~,~,~,~,J2]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,handles.din.fit_options);
            covar_x2=full(J2)\covariance_ratios/transpose(full(J2));%determine the covariance of numerically solved solution vector @ d2lam @ unique_n(dum)
            
%             cpd=@(x) [real(handles.din.solved.(name).complex_pred_freq-...
%                 master2(Dn_calc(d2lam_n2_calc(x(1),n1,unique_n(dum),x(2)),x(2)),Rliq_calc(delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,unique_n(dum))),unique_n(dum),x(3))).*sauerbrey(unique_n(dum),x(3)./1000));...
%                 imag(handles.din.solved.(name).complex_pred_freq-...
%                 master2(Dn_calc(d2lam_n2_calc(x(1),n1,unique_n(dum),x(2)),x(2)),Rliq_calc(delfnstar_liq_calc(zqliq_n2_calc(z_star_liq_n,dgliq_harm,unique_n(dum))),unique_n(dum),x(3))).*sauerbrey(unique_n(dum),x(3)./1000));...
%                 x(1)-x(1)];
%             [~,~,~,~,~,~,J_t]=lsqnonlin(cpd,solved,handles.din.lb,handles.din.ub,handles.din.fit_options);
%             J_t=full(J_t);
%             handles.din.solved.(name).J_cpd=J_t;%store the Jacobian for the complex predicted frequency
%             handles.din.stored_solutions.(name).error.J_cpd(count,:)={handles.din.cursor.selected_time,J_t};%store the Jacobian for the predicted frequency shift
%             
%             handles.din.stored_solutions.(name).error.cpf(count,:)=[handles.din.cursor.selected_time,...
%                 (2.*sqrt(J_t(1,1).^2.*covar_x2(1)+J_t(1,2).^2.*covar_x2(2,2)+J_t(1,3).^2.*covar_x2(3,3)+2.*J_t(1,1).*J_t(1,2).*covar_x2(1,2)+2.*J_t(1,1).*J_t(1,3).*covar_x2(1,3)+2.*J_t(1,2).*J_t(1,3).*covar_x2(2,3))),...
%                 (2.*sqrt(J_t(2,1).^2.*covar_x2(1)+J_t(2,2).^2.*covar_x2(2,2)+J_t(2,3).^2.*covar_x2(3,3)+2.*J_t(2,1).*J_t(2,2).*covar_x2(1,2)+2.*J_t(2,1).*J_t(2,3).*covar_x2(1,3)+2.*J_t(2,2).*J_t(2,3).*covar_x2(2,3)))];
            handles.din.stored_solutions.(name).error.d2lamn(count,:)=[handles.din.cursor.selected_time,sqrt(var_d2lam(unique_n(dum),n1,solved2(2),solved2(1),covar_x1)).*2];%d2lam error at unique_n(dum), twice the stdev
            handles.din.stored_solutions.(name).error.drho(count,:)=[handles.din.cursor.selected_time,sqrt(covar_x1(3,3)).*2.*1000];%drho error at unique_n(dum), twice the stdev (g/m^2)
            handles.din.stored_solutions.(name).error.grho(count,:)=[handles.din.cursor.selected_time,sqrt(var_grho_n1(unique_n(dum),solved2(2),solved2(1),solved(3),covar_x2)).*2./1000];%grho error at unique_n(dum), twice the stdev (Pa-g/cm^3)
            handles.din.stored_solutions.(name).error.phi(count,:)=[handles.din.cursor.selected_time,sqrt(covar_x1(2,2)).*2];%phi error, twice the stdev (deg.)
            handles.din.stored_solutions.(name).error.covar_x2(count,:)={handles.din.cursor.selected_time,{covar_x2}};
        end
    catch err_msg
        assignin('base','err_msg',err_msg);
        set(handles.status,'string','Status: Error in calculating errors for harm and diss ratios!','backgroundcolor','r','foregroundcolor','k');
        disp('!Error in calculating errors for harm and diss ratios!');
    end%try
    if get(handles.static_guess,'value')==0%only update the guess parameters for harmonic n1      
        name=['n_',num2str(unique_n(1)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
        guess_table=get(handles.uitable4,'data');
        guess_table{1,2}=handles.din.solved.(name).drho;
        guess_table{2,2}=handles.din.solved.(name).phi;
        guess_table{3,2}=handles.din.solved.(name).grho;
        guess_table{4,2}=handles.din.solved.(name).d2lam;
        set(handles.uitable4,'data',guess_table);
    end%get(handles.static_guess)==0
end%for dum=1:length(radiotot)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS MISC. FUNCTIONS
function [n1,n2,n3]=determine_harm_choice(harm_datasets)
%this function determines what harmonic datasets are being used to
%calculate the viscoelastic parameters.
%Example: 1,3:1 corresponds to n1=1, n2=3, n3=1
switch harm_datasets%determine which harmonics to use for the calculation
    case 1
        n1=1; n2=3; n3=1;
    case 2
        n1=1; n2=3; n3=3;
    case 3
        n1=1; n2=5; n3=1;
    case 4
        n1=1; n2=5; n3=5;
    case 5
        n1=3; n2=5; n3=3;
    case 6
        n1=3; n2=5; n3=5;
    case 7
        n1=1; n2=7; n3=7;
    case 8
        n1=3; n2=7; n3=7;
    case 9
        n1=5; n2=7; n3=7;
    case 10
        n1=1; n2=7; n3=1;
    case 11
        n1=3; n2=7; n3=3;
    case 12
        n1=5; n2=7; n3=5;
    case 13
        n1=1; n2=5; n3=3;
    case 14
        n1=5; n2=3; n3=3;
    case 15
        n1=3; n2=9; n3=3;
    case 16
        n1=5; n2=9; n3=5;
    case 17
        n1=5; n2=9; n3=9;
    case 18
        n1=3; n2=11; n3=11;
    case 19
        n1=3; n2=11; n3=3;
    case 20
        n1=7; n2=9; n3=7;
    case 21        
        n1=7; n2=9; n3=9;
    case 22
        n1=9; n2=11; n3=9;
    case 23
        n1=9; n2=11; n3=11;
    case 24
        n1=7; n2=11; n3=7;
    case 25
        n1=7; n2=11; n3=11;
    case 26
        n1=3; n2=5; n3=7;
    case 27
        n1=1; n2=3; n3=5;
    case 28
        n1=3; n2=5; n3=1;
end%switch harm_datasets

function [xstr,xdata]=determine_time_units(xdata,handles)
%This function determines what units the xaxis should be in and adjusts the
%xdata accordingly
%extract out the time array
    switch get(handles.time_units,'value')
        case 1%min
            xstr='Time (min.)';
        case 2%hr
            xstr='Time (hr)';
            xdata=xdata./60;%hr
        case 3%day
            xstr='Time (day)';
            xdata=xdata./(60.*24);%day
        case 4%Temperature
            xstr='Temperature (^oC)';
    end%switch get(handles.time_units,'value')


% --------------------------------------------------------------------
function pref_ClickedCallback(hObject, eventdata, handles)%under construction
if strcmp(hObject.State,'on')==1
    figure(99);clf(figure(99));
    p=uipanel('units','normalized','position',[0.05 0.55 0.9 0.45],'backgroundcolor','w','title','Solution contour plot options',...
        'fontweight','bold','fontsize',10);
    fit=uicontrol('style','pushbutton','units','normalized','string','Adjust fit options','position',[0.02 0.5 0.2 0.05]);
    set(fit,'callback',{@fit_options_fcn,handles});
else
    try delete(figure(99));end
end

function fit_options_fcn(~,~,handles)%under construction
% inspect(handles.din.fit_options);
    

% --- Executes on button press in dg_log.
function dg_log_Callback(hObject, eventdata, handles)
children=get(figure(3),'children');
if get(handles.dg_log,'value')==1
    set(children(2),'yscale','log');
else
    set(children(2),'yscale','linear'); 
end

% --- Executes on button press in error_values.
function error_values_Callback(hObject, eventdata, handles)
figure(99);clf(figure(99));pos=get(figure(99),'position');set(figure(99),'position',[pos(1)-100 pos(2) 1100 300],...
    'name','Covariance Matrix');
t=uitable('parent',figure(99),'units','normalized','position',[0.05 0.05 0.9 0.9]);
cnames=[{'del_f1'},{'del_g1'},{'del_f3'},{'del_g3'},{'del_f5'},{'del_g5'},{'del_f7'},{'del_g7'},{'del_f9'},{'del_g9'},{'del_f11'},{'del_g11'}];
rnames=[{'del_f1'},{'del_g1'},{'del_f3'},{'del_g3'},{'del_f5'},{'del_g5'},{'del_f7'},{'del_g7'},{'del_f9'},{'del_g9'},{'del_f11'},{'del_g11'}];
data=get(handles.error_values,'userdata');
set(t,'columnname',cnames,'rowname',rnames,'data',data);
set(t,'columneditable',[true true true true true true true true true true true true],'celleditcallback',{@store_error,handles});

function store_error(hObject,callbackdata,handles)
index=callbackdata.Indices;%extract which row and col the value was edited in
newdata=callbackdata.NewData;%extract the new data
set(handles.error_values,'userdata',get(hObject,'data'));


% --- Executes when entered data in editable cell(s) in uitable4.
function uitable4_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable4 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
try
    Index=eventdata.Indices;%extract index information that was edited
    if Index(1)==3
        handles.din.grho_state=0;
        handles.din.d2lam_state=1;
    elseif Index(1)==4
        handles.din.grho_state=1;
        handles.din.d2lam_state=0;
    else
    end%if Index(1)==3
    grho_state=handles.din.grho_state;
    d2lam_state=handles.din.d2lam_state;
    if grho_state==1&&d2lam_state==0
        handles=d2lam_guess_Callback(hObject, eventdata, handles);
    elseif grho_state==0&&d2lam_state==1
        handles=grho_guess_Callback(hObject, eventdata, handles);
    end%if grho_state==1&&d2lam_state==0
catch err_msg
    assignin('base','err_msg',err_msg);
    disp('Error in calculations!');
    set(handles.status,'string','Status: Error in calcualtions!','backgroundcolor','r','foregroundcolor','k');
end%try
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
handles.din.guess_label2=rand(1);
guidata(hObject,handles);

% --- Executes on button press in sim.
function sim_Callback(hObject, eventdata, handles)%This function adds a simulate button functionality
handles.din.simulate=1;
visco_table=get(handles.uitable2,'data');%extract the viscoelastic information
input_table=get(handles.uitable4,'data');%extract the input viscoelastic information
drho=input_table{1,2};%drho in g/m^2
phi=input_table{2,2};%phi in deg.
grho=input_table{3,2};%grho in Pa-g/cm^3
d2lam=input_table{4,2};%film thickness to shearwavelength ratio
f1=handles.din.constants.f1;%fundamentalr resonance freq in Hz
zq=handles.din.constants.zq;% the quartz load impedance
[radiotot,harmchoice]=radio_total(handles);%determine how many datasets and which harmonics that will be used to perform the calculation
for dum0=1:length(radiotot)
    handles.din.solved.harmchoice=harmchoice(dum0);
    %calculate the liquid load impedance
    dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
    dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance at dgliq_harm
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at
    [n1,n2,n3]=determine_harm_choice(harmchoice(dum0));%determine what harmonic datasets will be used to calculate viscoelastic parameters    
    unique_n=1:2:11;
    for dum=1:length(unique_n)
        name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
        handles.din.solved.(name).drho=drho;
    end
    for dum=1:length(unique_n)        
        handles=calc_viscoelastic_par(handles,unique_n(dum),f1,z_star_liq_n,zq,[d2lam phi drho],dgliq_harm,name,n1);%calculate the viscoelastic parameters
        %store the calculated values in a matrix that will be used to output the results in the handles.uitable2
        calc_table((unique_n(dum)+1)/2,:)=[handles.din.solved.(name).drho,handles.din.solved.(name).grho,handles.din.solved.(name).phi,...
            handles.din.solved.(name).d2lam,handles.din.solved.(name).lam_rho,real(handles.din.solved.(name).norm_delfstar),...
            imag(handles.din.solved.(name).norm_delfstar),abs(handles.din.solved.(name).cf),handles.din.solved.(name).del2d];
        %store the predicted freq shifts in a matrix that will be used to output the results in the handles.uitable3
        pred_shifts((unique_n(dum)+1)/2,:)=[nan,real(handles.din.solved.(name).complex_pred_freq),...
            nan,imag(handles.din.solved.(name).complex_pred_freq)];
    end
    %update the data tables int he gui
    set(handles.uitable2,'data',calc_table);%update the viscoelastic parameter table
    set(handles.uitable3,'data',pred_shifts);%update the predicted frequency shift table     
    set(handles.status,'string',['Status: Simulation complete for ',num2str(n1),num2str(n2),num2str(n3)]);
end
handles.din.simulate=0;
guidata(handles.figure1,handles);


% --- Executes on button press in fit_all.
function fit_all_Callback(hObject, eventdata, handles)
if isfield(handles.raw,'filename')==1&&isfield(handles.raw,'pathname');
    disp(['Raw spectras extracted from : ', handles.raw.pathname,handles.raw.filename]);
    initial_time=str2double(handles.initial_time.String);
    solve_inc=str2double(handles.solve_inc.String);
    extent=str2double(handles.extent.String);
    harm=floor(str2double(inputdlg('Which harmonic do you wish to refit (integers only!)','Choose harmonic to refit',1,{'3'})));
    harm_dataf=handles.din.(['harmf',num2str(harm)]);
%     harm_datag=handles.din.(['harmg',num2str(harm)]);
    if strcmp(class(harm),'double')==0||isempty(harm)==1||isnan(harm)==1
        disp('Harmonic input was not a number, not odd, and/or less than 12!');
        return
    else
        disp(['Harmonic to be fitted: ',num2str(harm)]);
    end
    disp(['Fitting from global index [',num2str(initial_time),'] to [',num2str(extent),'] in increments of [',num2str(solve_inc),']']);    
        for dum=initial_time:solve_inc:extent
            try
                handles=guidata(handles.figure1);
                timeline=handles.din.rawdata.freq_shift(:,1);
                timepoint0=timeline(dum)
                event.Position(1)=timepoint0;
                handles=rawfit(handles,harm,event);
                drawnow;
                refit_panel=findall(0,'type','figure','name','Refit panel');
                h=guidata(refit_panel);
                try
                    h.guess_values_options.Value=4;
                catch
                    rpc=refit_panel.Children;
                    h.guess_values_options=findall(rpc,'style','popupmenu');
                    h.guess_values_options.Value=4;
                end
                handles.raw.index=find(harm_dataf(:,1)==timepoint0);
                guidata(handles.figure1,handles);            
                for repeat=1:1
                    fit_button_callback(hObject,1,guidata(handles.figure1),harm,h.guess_values_options,h.a4,h.a1,h.a2,h.a3);
                end
                drawnow;
                handles=guidata(handles.figure1);    
                guidata(handles.figure1,handles);            
                accept_fcn(1,1,handles);
                drawnow;
            catch err
                assignin('base','err',err);
                disp(['Error! Index: ',num2str(dum)]);        
            end
        end
        disp('Automated refitting process finished.');    
else
    disp('Raw spectras have not been loaded!');
end


% --------------------------------------------------------------------
function txt=myupdatefcn(~,event,handles)
handles=guidata(handles.figure1);
dis=get(event.Target,'displayname');
if length(dis)==3
    harm=str2double(dis(end));
elseif length(dis)==4
    harm=str2double(dis(end-1:end));
end
if strcmp(get(get(event.Target,'parent'),'tag'),'axes1')
    str='\Delta f/n: ';
elseif strcmp(get(get(event.Target,'parent'),'tag'),'axes2')
    str='/Delta/Gamma: ';
end
dis_t=get(event.Target,'xdata');
index=find(dis_t==event.Position(1));
index2=find(handles.din.rawdata.freq_shift(:,1)==event.Position(1));
txt={['Time: ',num2str(event.Position(1)),' min'],...
    [str,num2str(event.Position(2)),' Hz'],get(event.Target,'displayname'),...
    ['Local Index: ',num2str(index)],...
    ['Global Index: ',num2str(index2)]};
handles.raw.index=index;
guidata(handles.figure1,handles);

% --------------------------------------------------------------------
function import_cov_ClickedCallback(hObject, eventdata, handles)
try
    [filename,filepath,~]=uigetfile('.mat','Load covariance matrix',handles.din.filepath);%get the filepath and filename
    if isstr(filename)==0
        return
    end
    set(handles.status,'string','Status: Importing covariance matrix...','backgroundcolor','k','foregroundcolor','r');
    covar_data=load([filepath,filename]);
    handles.error_values.UserData=covar_data.covar1;
    guidata(handles.figure1,handles);
    disp(['Status: Covariance matrix imported! filename: ',[filepath,filename]]);
    set(handles.status,'string',['Status: Covariance matrix imported! filename: ',[filepath,filename]],'backgroundcolor','k','foregroundcolor','r');    
    set(handles.error_values,'tooltipstring',['filename: ',filepath,filename]);
catch err_msg
    assignin('base','err_msg','err_msg');
    disp('Import failed!');
    set(handles.status,'string','Status: Import failed!','backgroundcolor','r','foregroundcolor','k');
end
