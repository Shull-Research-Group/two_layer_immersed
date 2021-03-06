function varargout = two_layer_immersed_v003(varargin)
% Shull research group. Original author: Chyi-Huey Joshua Yeh
%This version is currently in beta mode
% TWO_LAYER_IMMERSED_V003 MATLAB code for two_layer_immersed_v003.fig
%      TWO_LAYER_IMMERSED_V003, by itself, creates a new TWO_LAYER_IMMERSED_V003 or raises the existing
%      singleton*.
%
%      H = TWO_LAYER_IMMERSED_V003 returns the handle to a new TWO_LAYER_IMMERSED_V003 or the handle to
%      the existing singleton*.
%
%      TWO_LAYER_IMMERSED_V003('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWO_LAYER_IMMERSED_V003.M with the given input arguments.
%
%      TWO_LAYER_IMMERSED_V003('Property','Value',...) creates a new TWO_LAYER_IMMERSED_V003 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before two_layer_immersed_v003_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to two_layer_immersed_v003_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help two_layer_immersed_v003

% Last Modified by GUIDE v2.5 13-Jul-2016 13:07:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @two_layer_immersed_v003_OpeningFcn, ...
                   'gui_OutputFcn',  @two_layer_immersed_v003_OutputFcn, ...
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


% --- Executes just before two_layer_immersed_v003 is made visible.
function two_layer_immersed_v003_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to two_layer_immersed_v003 (see VARARGIN)

% Choose default command line output for two_layer_immersed_v003
handles.output = hObject;

%set default values or settings
handles=home_state(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes two_layer_immersed_v003 wait for user response (see UIRESUME)
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
handles.din.fit_options=optimset('Display','off','tolfun',1e-11,'tolx',1e-11,'Maxiter',10000000,'MaxFunEvals',10000000);
set(handles.viewraw,'userdata',[],'tooltipstring',['view raw spectra']);
set(handles.harm1,'value',1);
set(handles.harm3,'value',1);
set(handles.harm5,'value',1);
set(handles.radio2,'value',1);
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
dcm=datacursormode(handles.figure1);
set(dcm,'UpdateFcn',{@myupdatefcn,handles,handles.viewraw},'enable','off');
for dum=1:1:6
    name=['harmchoice',num2str(dum)];
    set(handles.(name),'string',[{'1,3:1'},{'1,3:3'},{'1,5:1'},{'1,5:5'},{'3,5:3'},{'3,5:5'},...
        {'1,7:7'},{'3,7:7'},{'5,7:7'},{'1,7:1'},{'3,7:3'},{'5,7:5'},{'1:5,3'},{'5,3:3'},{'3,9:3'},...
        {'5,9:5'},{'5,9:9'},{'3,11:11'},{'3,11:3'},{'7,9:7'},{'7,9:9'},{'9,11:9'},{'9,11:11'},...
        {'7,11:7'},{'7,11:11'},{'3,5:7'},{'1,3:5'}]);
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
function varargout = two_layer_immersed_v003_OutputFcn(hObject, eventdata, handles) 
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
    handles=import_data(handles,filename,filepath);
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
if strcmp(handles.viewraw.State,'off')
    handles.din.extent=length(timepoints);
    set(handles.extent,'string',handles.din.extent);
end
initial_time_Callback(handles.initial_time, 1, handles);
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
res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
phi_var=linspace(0,90,res);%define range of phi vales
d2lam_var=linspace(0,0.2,res);%define range of d2lam values
harm_ratio_calc=handles.din.contour.harm_ratio_calc;
diss_ratio_calc=handles.din.contour.diss_ratio_calc;
f1=handles.din.constants.f1;%fundamental resonance freq
zq=handles.din.constants.zq;%quartz  load impedance
%determine what n1, n2, n3 should be based on the value of the
%handles.contour_choice handle
[n1,n2,n3]=determine_harm_choice(get(handles.contour_choice,'value'));
%calculate experimental harmonic and dissipation ratio
del_f_n2=handles.din.cursor.(['interp_harmfi',num2str(n2)]);%Df_n2
del_f_n1=handles.din.cursor.(['interp_harmfi',num2str(n1)]);%Df_n1
del_g_n3=handles.din.cursor.(['interp_harmgi',num2str(n3)]);%Dg_n3
del_f_n3=handles.din.cursor.(['interp_harmfi',num2str(n3)]);%Df_n3
harm_ratio_exp=(del_f_n2.*n1)./(del_f_n1.*n2);%experimental harmonic ratio
diss_ratio_exp=del_g_n3./del_f_n3;%experimental dissipation ratio
%calculate the error in the harm and diss ratios
[harm_ratio_error,diss_ratio_error,handles]=calc_ratio_error(handles,n1,n2,n3,del_f_n2,del_f_n1,del_g_n3,del_f_n3);
%calculate resonance frequencies at n1, n2, and n3
f_n1=f1*n1;%resonant frequency at n1 in Hz
f_n2=f1*n2;%resonant frequency at n2 in Hz
f_n3=f1*n3;%resonant frequency at n3 in Hz
guess_table=get(handles.uitable4,'data');
if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled    
    try
        drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
        drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
        drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
        drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values
        guess_table{1,2}=drho_ave;
        set(handles.uitable3,'data',guess_table);
    catch
        drho_ave=guess_table{1,2}./1000;%kg/m^2
    end
else%if the manual guess for drho has been enabled
%     drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
    drho_ave=guess_table{1,2}./1000;%kg/m^2
end%if get(handles.edit_drho,'value')==0
dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
n=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at  
z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance
z_star_liq_n1=zqliq_n2_calc(z_star_liq_n,n,n1);%liq load impedance at harmonic n1
z_star_liq_n2=zqliq_n2_calc(z_star_liq_n,n,n2);%liq load impedance at harmonic n2
z_star_liq_n3=zqliq_n2_calc(z_star_liq_n,n,n3);%liq load impedance at harmonic n3
if handles.din.contour.label~=handles.din.guess_label%only regenerate the ratio matrices if the parameters have changed
    h=waitbar(0,'Generating plot...');%create progress bar
    for dum1=1:length(phi_var);
        try
            waitbar(dum1/length(phi_var),h);
        catch
            set(handles.status,'string','Cancelled plot generation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled plot generation.');
            return
        end%try
        for dum2=1:length(d2lam_var)       
            if get(handles.d2lam_choice,'value')==2 %relative to d2lambda_n2
                d2lam_n2=d2lam_var(dum2);%d2lambda at harmonic n2
                d2lam_n1=d2lam_n2_calc(d2lam_n2,n2,n1,phi_var(dum1));% calc. d2lambda ratio at harmonic n1
                d2lam_n3=d2lam_n2_calc(d2lam_n1,n1,n3,phi_var(dum1));% calc. d2lambda at harmonic n3            
                handles.din.n_ref=n2;
            elseif get(handles.d2lam_choice,'value')==1 %relative to d2lambda_n1
                d2lam_n1=d2lam_var(dum2);% d2lambda at harmonic n1
                d2lam_n2=d2lam_n2_calc(d2lam_n1,n1,n2,phi_var(dum1));% calc. d2lambada ratio at harmonic n2
                d2lam_n3=d2lam_n2_calc(d2lam_n1,n1,n3,phi_var(dum1));% calc. d2lambda at harmonic n3
                handles.din.n_ref=n1;
            end% if get(handles.d2lam_choice,'value')==2
            %only used drho_n2 because drho_n1, n2, n3 should be similar
            harm_ratio_calc(dum1,dum2)=real(master(d2lam_n2,phi_var(dum1),z_star_liq_n2,f_n2,drho_ave))./real(master(d2lam_n1,phi_var(dum1),z_star_liq_n1,f_n1,drho_ave));%calculate the harmonic ratio
            diss_ratio_calc(dum1,dum2)=imag(master(d2lam_n3,phi_var(dum1),z_star_liq_n3,f_n3,drho_ave))./real(master(d2lam_n3,phi_var(dum1),z_star_liq_n3,f_n3,drho_ave));%calculate the dissipation ratio
        end%for dum2=length(d2lam_var)
    end%for dum1=length(phi_var);
end%if handles.contour.label~=handles.din.guess_label

%store the calculated ratio in the handles structure so that the plot does
%not have to be reproduced every single time the callback function runs.
handles.din.contour.harm_ratio_calc=harm_ratio_calc;%store the harm_ratio_calc
handles.din.contour.diss_ratio_calc=diss_ratio_calc;%store the diss_ratio_calc
handles.din.contour.label=rand(1);%store a "label" for the stored calc. ratios

%calculate the solutions
if get(handles.qcm_map_sol,'value')==1
    %find indices that fall inbetween the lower and upperbound of the
    %experimentally determined ratios
    %pre-allocate variables
    d2lam_hsol=[];
    phi_hsol=[];
    d2lam_dsol=[];
    phi_dsol=[];
    d2lam_dhsol=[];
    phi_dhsol=[];
    try
        waitbar(0,h,'calculating harm ratio solution...');
    catch
        h=waitbar(0,'calculating harm ratio solution...');
    end%try
    %find row and col indices that corresponds to ratio values that fall within error of the experimentally determined ratio values
    [rh,ch]=find(harm_ratio_calc>=harm_ratio_exp-harm_ratio_error&harm_ratio_calc<=harm_ratio_exp+harm_ratio_error);%harm ratio
    [dr,dc]=find(diss_ratio_calc>=diss_ratio_exp-diss_ratio_error&diss_ratio_calc<=diss_ratio_exp+diss_ratio_error);%diss ratio
    [rdh,cdh]=find(harm_ratio_calc>=harm_ratio_exp-harm_ratio_error&harm_ratio_calc<=harm_ratio_exp+harm_ratio_error&...
        diss_ratio_calc>=diss_ratio_exp-diss_ratio_error&diss_ratio_calc<=diss_ratio_exp+diss_ratio_error);%harm and diss ratio
    %find the coresponding d2lam and phi values that will give the
    %target harm ratio
    for dum1=1:length(rh)%harm ratio
        try%update the status bar
            waitbar(dum1./length(rh),h);
        catch
            set(handles.status,'string','Cancelled solution calculation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled solution calculation.');
            return
        end%try
        d2lam_hsol(dum1)=d2lam_var(ch(dum1));%find corresponding d2lam values that will give the target harm ratio
        phi_hsol(dum1)=phi_var(rh(dum1));%find corresponding phi values that will give the target harm ratio
    end%for dum=1:length(rh)    
    try%update the status bar
        waitbar(0,h,'calculating diss ratio solution...');
    end%try
    %find the corresponding d2lam and phi values that will give the target
    %diss ratio
    for dum1=1:length(dr)%diss ratio
        try%update the status bar
            waitbar(dum1./length(dr),h);
        catch
            set(handles.status,'string','Cancelled solution calculation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled solution calculation.');
            return
        end%try
        d2lam_dsol(dum1)=d2lam_var(dc(dum1));%find corresponding d2lam values that will give the target diss ratio
        phi_dsol(dum1)=phi_var(dr(dum1));%find corresponding phi values that will give the target harm ratio
    end%for dum=1:length(rh)    
    try%update the status bar
        waitbar(0,h,'calculating diss and harm ratio solution...');
    end%
    %find the corresponding d2lam and phi values that will give the target
    %harm and diss ratio
    for dum1=1:length(rdh)%diss ratio and harm ratio
        try%update the status bar
            waitbar(dum1./length(rdh),h);
        catch
            set(handles.status,'string','Cancelled solution calculation. Ready...','backgroundcolor','k','foregroundcolor','r');
            disp('Cancelled solution calculation.');
            return
        end%try
        d2lam_dhsol(dum1)=[d2lam_var(cdh(dum1))];%find corresponding d2lam values that will give the target diss and harm ratio
        phi_dhsol(dum1)=[phi_var(rdh(dum1))];%find corresponding phi values that will give the target diss and harm ratio
    end%for dum=1:length(rh)  
    try  delete(h); end%try to delete the status bar
end%if get(handles.qcm_map_sol,'value')==1

if get(handles.qcm_map_sol,'value')==1
    handles.din.contour.d2lam_hsol=d2lam_hsol;%store the d2lam values that will give the experimental harm ratio
    handles.din.contour.phi_hsol=phi_hsol;%store the phi values that will give the experimental harm ratio
    handles.din.contour.d2lam_dsol=d2lam_dsol;%store the d2lam values that will give the experimental diss ratio
    handles.din.contour.phi_dsol=phi_dsol;%store the phi values that will give the experimental diss ratio
    handles.din.contour.d2lam_dhsol=d2lam_dhsol;%store the d2lam values that will give the experimental harm and diss ratio
    handles.din.contour.phi_dhsol=phi_dhsol;%store the phi values that will give the experimental harm and diss ratio
else
    try   delete(h); end%try to delete the staus bar
end%if get(handles.qcm_map_sol,'value')==1
handles.din.guess_label=handles.din.contour.label;%rewrite the label for the guess parameters
guidata(hObject,handles);

%plot the contour map
ff=figure();
pos=get(ff,'position');
if pos(3)<=800
    set(ff,'position',[pos(1)*.5 pos(2) pos(3)*2 pos(4)]);
end%if pos(3)<=0.71
drawnow;
s1=subplot(1,2,1);%plot the harmonic ratio map
contourf(d2lam_var,phi_var,harm_ratio_calc,linspace(0,1,res),'edgecolor','none');
xlabel(s1,['d/\lambda_',num2str(handles.din.n_ref)],'fontweight','bold','fontsize',14);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',14);
t1=title([num2str(n1),'\Deltaf_',num2str(n2),'/',num2str(n2),'\Deltaf_',num2str(n1)],'fontweight','bold','fontsize',14);
c1=colorbar;
set(c1,'fontsize',12);
set(s1,'fontsize',14);
drawnow;
s2=subplot(1,2,2);%plot the dissipation ratio map
contourf(d2lam_var,phi_var,diss_ratio_calc,linspace(-1.1,0,res),'edgecolor','none');
xlabel(s2,['d/\lambda_',num2str(handles.din.n_ref)],'fontweight','bold','fontsize',14);
ylabel(s2,'\phi (deg.)','fontweight','bold','fontsize',14);
t2=title(['\Delta\Gamma_',num2str(n3),'/\Deltaf_',num2str(n3)],'fontweight','bold','fontsize',14);
c2=colorbar;
set(c2,'fontsize',12);
set(s2,'fontsize',14);
linkaxes([s1,s2],'xy');%link the axes
drawnow;
set(handles.status,'string','Status: Plot succesfully generated!','foregroundcolor','r','backgroundcolor','k');

%plot the solutions
if get(handles.qcm_map_sol,'value')==1
hold(s1,'on');
contour(s1,d2lam_var,phi_var,harm_ratio_calc,[harm_ratio_exp-harm_ratio_error harm_ratio_exp+harm_ratio_error],'edgecolor','k','linestyle','-');
contour(s1,d2lam_var,phi_var,diss_ratio_calc,[diss_ratio_exp-diss_ratio_error diss_ratio_exp+diss_ratio_error],'edgecolor','k','linestyle','--');
hold(s2,'on');
contour(s2,d2lam_var,phi_var,harm_ratio_calc,[harm_ratio_exp-harm_ratio_error harm_ratio_exp+harm_ratio_error],'edgecolor','k','linestyle','-');
contour(s2,d2lam_var,phi_var,diss_ratio_calc,[diss_ratio_exp-diss_ratio_error diss_ratio_exp+diss_ratio_error],'edgecolor','k','linestyle','--');
t1.String=[t1.String,' = ',num2str(harm_ratio_exp),'\pm',num2str(harm_ratio_error)];
t2.String=[t2.String,' = ',num2str(diss_ratio_exp),'\pm',num2str(diss_ratio_error)];
end%if get(handles.qcm_map_sol,'value')==1
linkaxes([s1 s2],'xy');
drawnow
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
ff1=figure(4);clf(ff1);pos=get(figure(4),'position');
set(ff1,'position',[pos(1),pos(2),1466,414]);
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
ff2=figure(3);clf(ff2);pos=get(figure(3),'position');
set(ff2,'position',[pos(1),pos(2),1466,414]);
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


function qcm_map_Callback(hObject, eventdata, handles)
%this function plots out the contour plots for the master equation
set(handles.status,'string','Status: Generating QCM maps...','backgroundcolor','k','foregroundcolor','r');
disp('Generating QCM maps...');
qcm_map=handles.din.qcm_map;
res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
phi_var=linspace(0,90,res);%define range of phi vales
d2lam_var=linspace(0,0.01,res);%define range of d2lam values
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
ff=figure;
pos=get(ff,'position');
if pos(3)<=800
    set(ff,'position',[pos(1)*.5 pos(2) pos(3)*2 pos(4)]);
end%if pos(3)<=0.71
s1=subplot(1,2,1);%plot the real component
contourf(d2lam_var,phi_var,real(qcm_map),linspace(-1.1,-.9,300),'edgecolor','none');hold on;
xlabel(s1,'d/\lambda_n','fontweight','bold','fontsize',12);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Deltaf_n','/\Deltaf_s_',num2str(n)],'fontweight','bold','fontsize',12);
colormap(jet); colorbar;
waitbar(.95,h);
s2=subplot(1,2,2);%plot the imaginary component
contourf(d2lam_var,phi_var,imag(qcm_map),linspace(0,0.1,300),'edgecolor','none');hold on;
xlabel(s2,'d/\lambda_n','fontweight','bold','fontsize',12);
ylabel(s2,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Delta\Gamma_n','/\Deltaf_s_',num2str(n)],'fontweight','bold','fontsize',12);
colormap(jet); colorbar;
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
    edata1=abs([ydata1-handles.din.stored_solutions.(temp).error.d2lam(:,2),handles.din.stored_solutions.(temp).error.d2lam(:,3)-ydata1]);%calculate the error for drho
    edata2=abs([ydata2-handles.din.stored_solutions.(temp).error.phi(:,2),handles.din.stored_solutions.(temp).error.phi(:,3)-ydata2]);%calculate the error for phi
    %need to make sure that the arrays have the same dimensions
    ind1=find(edata2(:,1)~=0);
    try ydata1=ydata1(ind1);end
    try ydata2=ydata2(ind1);end
    try edata1=edata1(ind1,:);end
    try edata2=edata2(ind1,:);end
    k=strfind(temp,'_');
    L_string=[L_string,{temp(k(end)+1:end)}];%string for legend box
    try
        errorbar(s1,ydata1,ydata2,edata1(:,1),edata1(:,2),...
            'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        errorbar(s2,ydata1,ydata2,edata2(:,1),edata2(:,2),...
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
    edata1=abs([ydata1-handles.din.stored_solutions.(temp).error.drho(:,2),handles.din.stored_solutions.(temp).error.drho(:,3)-ydata1]);%calculate the error for drho
    edata2=abs([ydata2-handles.din.stored_solutions.(temp).error.grho(:,2),handles.din.stored_solutions.(temp).error.grho(:,3)-ydata2]);%calculate the error for grho
    edata3=abs([ydata3-handles.din.stored_solutions.(temp).error.phi(:,2),handles.din.stored_solutions.(temp).error.phi(:,3)-ydata3]);%calculate the error for phi
    %need to make sure that the arrays have the same dimensions
    ind1=find(edata2(:,1)~=0);
    xdata=xdata(ind1);
    try ydata1=ydata1(ind1);end
    try ydata2=ydata2(ind1);end
    try ydata3=ydata3(ind1);end
    try edata1=edata1(ind1,:);end
    try edata2=edata2(ind1,:);end
    try edata3=edata3(ind1,:);end
    if get(handles.norm_time,'value')==1%run this code if the option to normalize the timepoints is turned on
        xdata=xdata-min(xdata);
    end%if get(handles.norm_time,'value')==1
    k=strfind(temp,'_');
    L_string=[L_string,{temp(k(end)+1:end)}];%string for legend box
    try
        errorbar(s1,xdata,ydata1,edata1(:,1),edata1(:,2),'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        errorbar(s2,xdata,ydata2,edata2(:,1),edata2(:,2),'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
        errorbar(s3,xdata,ydata3,edata3(:,1),edata3(:,2),'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s3,'on');
    catch
        set(handles.status,'string',['Status: Could not plot errorbar dataset for ',temp(end-2:end)],'backgroundcolor','y','foregroundcolor','r');
        disp(['Could not plot errorbar dataset for ',temp(end-2:end)]);
        plot(s1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s1,'on');
        plot(s2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s2,'on');
        plot(s3,xdata,ydata3,'linewidth',2,'markersize',14,'marker',handles.din.marker_style{dum},'color',handles.din.marker_color{dum},'linestyle','none');hold(s3,'on');
    end%try        
end%for dum=1:length(index2)


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
        ldata1=ydata1-real(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,2));
        udata1=real(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,3))-ydata1;
        ydata2=handles.din.stored_solutions.(temp).g_shifts(:,3);%Calc Dg       
        ldata2=ydata2-imag(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,2));
        udata2=imag(handles.din.stored_solutions.(temp).error.complex_pred_freq(:,3))-ydata2;
%         J=find(ldata2<0)
%         ldata2(J)=1e-2;
%         j=find(udata2<0)
%         udata2
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
%             errorbar(r1,xdata(1),ydata1(1),ldata1(1),udata1(1),'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none');
%             errorbar(r2,xdata(1),ydata2(1),ldata2(1),udata2(1),'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none'); 
        elseif flag==1
            index3=0.5*(dum1+1);
            hold(r1,'on');hold(r1,'on');
%             plot(r1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
%             plot(r2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
            errorbar(r1,xdata,ydata1,ldata1,udata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
            errorbar(r2,xdata,ydata2,ldata2,udata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none'); 
        end%if flag==0
    end%for dum=1:length(index2)
end%for dum1=1:2:11



function qcm_map_sol_Callback(hObject, eventdata, handles)



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
if strcmp(get(handles.viewraw,'enable'),'off')
disp('loading harmonic...');
end
set(handles.status,'string','Loading harmonic...');
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
if get(hObject,'value')==1&&strcmp(get(handles.viewraw,'enable'),'off')
    h=waitbar(0,'Loading harmonic...');
    for dum=1:3;
        pause(1);waitbar(dum/3,h);
    end;delete(h);
end
if strcmp(get(handles.viewraw,'enable'),'off')
disp('Harmonic loaded');
end
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
    set(hObject,'tooltipstring',['Initial timepoint: ',num2str(timepoints(ind1))]);
else
    set(hObject,'tooltipstring',['Initial timepoint: ',num2str(timepoints(ind1)-min(timepoints))]);
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
    set(hObject,'tooltipstring',['Initial timepoint: ',num2str(timepoints(ind1))]);
else
    set(hObject,'tooltipstring',['Initial timepoint: ',num2str(timepoints(ind1)-min(timepoints))]);
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
    end%dum2=1:length(handles.din.harmtot)    
    handles=qcm_solver(handles,count);%solve for viscoelastic parameters  
    count=count+1;
end%for dum=length(timepoints)
try delete(h); end%try to delete the status bar
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

function d2lam_n2=d2lam_n2_calc(d2lam_n1,n1,n2,phi)
%this function calculates the d2lambda ratio for harmonic n2 and puts it in
%the d2lam_n2 variable. The inputs are as follows: 
%d2lam_n1: d2lambda ratio at harmonic n1
%2) n1: harmonic n1
%3) n2: harmonic n2
%4) phi: viscoelastic phase angle
d2lam_n2=(d2lam_n1).*(n1./n2).^((phi./180)-1);

function grho=grho_calc(drho,d2lam,f_n,phi)
grho=(drho.*f_n.*cosd(phi./2)./d2lam).^2;%in Pa-kg/m^3

function z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq)
%This function calculates the liquid load impedance based on the small load
%approximation impedance equation at harmonic n. The inputs are as follows:
%dfliq_n: the frequency shift on the location of the resonance peak relative to a
%bare xtal in air at harmonic n
%dgliq_n: the dissipation shift relative to a bare xtal in air at harmonic n
%f1: the fundamental resonance frequency
%zq: the load impedance of quartz
z_star_liq_n=((dfliq_n+1i.*dgliq_n).*pi.*zq)./(1i.*f1);

function z_star_liq_n2=zqliq_n2_calc(zqliq_n1,n1,n2)
%This function calculates the liquid load impedance at harmonic n2 based on
%the Kanazawa and Gordon relation. The inputs are as follows:
%zqliq_n1: zqliq_n1
%n1: harmonic n1
%n2: harmonic n2
z_star_liq_n2=zqliq_n1.*sqrt(n2./n1);

function drho=drho_est(df,n,f1,zq)
%This funtion estimates the drho value based on the Sauerbrey equation. The
%inputs are as follows:
%df:  frequency shift at the nth harmonic
%n: harmonic of the freq shift
%f1: fundametnal resonant frequency
%zq: load impedance of quartz
drho=-(df.*zq)./(2.*n.*f1.^2);

function drho_recalc=drho_calc(df,zq,f1,n,norm_delfstar)
%This function recalculates the areal mass, drho, with the parameters that
%were calculated using the estimated drho value.
drho_recalc=-drho_est(df,n,f1,zq)./real(norm_delfstar);

function norm_delfstar=master(d2lam_n,phi_n,z_star_liq_n,f_n,drho_n)
%This function represents the complex frequecy shift normalized by the
%Sauerbrey equation. Thei inputs are as follows:
%d2lam_n: the s2lambda ratio at harmonic n
%phi_n: the viscoelastic phase angle at harmonic n
%z_star_liq_n: the load impedance of water at harmonic n
%f_n: the resonant frequncy at harmonic n
%drho: the areal mass
norm_delfstar=-(((1./(d2lam_n.*(1-1i.*tand(phi_n./2)))).^2)-((z_star_liq_n./(f_n.*drho_n)).^2))./...
    (2.*pi.*(((cot(2.*pi.*d2lam_n.*(1-1i.*tand(phi_n./2))))./(d2lam_n.*(1-1i.*tand(phi_n./2))))+((1i.*z_star_liq_n)./(drho_n.*f_n))));

function fsn=sauerbrey(n,f1,drho,zq)
%This function calculates the sauerbbrey frequency shift
%Note that fsn notation is defined by Denolf et al., Langmuir, 2011 paper
fsn=(2.*n.*(f1.^2).*drho)./(zq);

function etarho=etarho_calc(dgliq_n,zq,n,f1)
%This function calculates eta-rho or the density times the viscosity of the
%fluid.
%dgliq_n: The dissipation shift of immersed bare QCM crystal in the liquid
%n: harmonic order
etarho=((dgliq_n.^2).*(zq.^2).*pi)./(n.*(f1.^3));

function delfnstar_liq=delfnstar_liq_calc(f1,zq,z_star_liq_n)
%this function calculates the complex frequency shift of a bare qcm xtal in
%liquid
%Note that delfnstar_liq is still dependent on n since z_star_liq_n is
%dependent on n.
delfnstar_liq=(1i.*f1.*z_star_liq_n)./(pi.*zq);


function Rliq=Rliq_calc(delfnstar_liq,n,drho,f1,zq)
%dfliq_n: freq shift of bare xtal in liq
%dgliq_n: diss shift of bare xtal in liq
%n: harmonic order
%drho: areal mass
Rliq=(delfnstar_liq)./sauerbrey(n,f1,drho,zq);

function Dn=Dn_calc(d2lam,phi)
%n: harmonic order
%d2lam_n: d/lambda
%phi: phase angle in degrees
Dn=2.*pi.*d2lam.*(1-1i.*tand(phi./2));

function delfstar2layer=master2(Dn, Rliq)
%Dn: see function Dn
%Rliq: see function Rliq
delfstar2layer=-((Dn.^-2)+(Rliq^2))./((cot(Dn)/Dn)+(Rliq));

function harm_ratio_error=h_ratio_error(s_fn1,s_fn2,cov_fn1_fn2,del_f_n1,del_f_n2,n1,n2)
%this function calculates the error assosciated with the harmonic ratio.
%This error calculation is based on standard propagation of error. See
%DeNolf et al. Langmuir, 2014 for more details on error analysis.
harm_ratio_error=sqrt(((n2./n1).^2).*(1./(del_f_n2.^2)).*((s_fn1.^2)-2.*((del_f_n1./del_f_n2).^2).*cov_fn1_fn2+((del_f_n1./del_f_n2).^2).*(s_fn2.^2)));

function diss_ratio_error=d_ratio_error(s_gn3,s_fn3,cov_gn3_fn3,del_g_n3,del_f_n3)
%This function calcualtes the error associated with the dissipation ratio.
%This error calculation is based on standard propgation of error. See
%DeNolf et al. Langmuir, 2014 for more details on error analysis.
diss_ratio_error=(1./(del_f_n3.^2)).*(((del_g_n3./del_f_n3).^2).*(s_fn3.^2)-2.*((del_g_n3./del_f_n3).^2).*cov_gn3_fn3+((s_gn3).^2));

function [harm_ratio_error,diss_ratio_error,handles]=calc_ratio_error(handles,n1,n2,n3,del_f_n2,del_f_n1,del_g_n3,del_f_n3)
%This function calculates the error for the harmonic and dissipation error
%based in the experimental error from the frequency shifts.
%handles: gui handles structure
%n1: first harmonic dataset used to calculate the harmonic ratio
%n2: second harmonic dataset used to calculate the harmonic ratio
%n3:  third harmonic dataset used to calcualte the dissipation ratio
%calculate error in ratios
data=get(handles.error_values,'userdata');
errors=diag(data);
s_fn1=errors(n1);
s_fn2=errors(n2);
s_fn3=errors(n3);
s_gn1=errors(n1+1);
s_gn2=errors(n2+1);
s_gn3=errors(n3+1);
cov_fn1_fn2=data(n1,n2);
cov_gn3_fn3=data(n3+1,n3);
harm_ratio_error=h_ratio_error(s_fn1,s_fn2,cov_fn1_fn2,del_f_n1,del_f_n2,n1,n2);%error in the harmonic ratio
diss_ratio_error=d_ratio_error(s_gn3,s_fn3,cov_gn3_fn3,del_g_n3,del_f_n3);%error in the dissipation ratio


function F=immersed_2layer_solver(handles,n1,n2,n3,n,z_star_liq_n,harm_ratio_exp,diss_ratio_exp,guess_values)
%This function is the function that is called by fsolve. It calculates the
%dissipation and harmonic ratio and takes the difference relative to the
%experimental values (F).
%n1: harmonic at n1
%n2: harmonic at n2
%n3: harmonic at n3
%n: hamonic at n or the harmonic in which the liq load impedance is
%associated with
%z_star_liq_n: liq load impedance at harmonic n
%harm_ratio_exp: experimentally obtained harmonic ratio
%diss_ratio_exp: experimentally obtained dissipation ratio
%guess_values(1): d2lam_n1 guess
%guess_values(2): phi
%guess_values(3): drho
f1=handles.din.constants.f1;
zq=handles.din.constants.zq;
d2lam_n1=guess_values(1);%guess d2lambda ratio at harmonic n1
phi=guess_values(2);%guess viscoelastic phase angle
drho=guess_values(3);%guess drho value
d2lam_n2=d2lam_n2_calc(d2lam_n1,n1,n2,phi);%based on guess, calc. d2lambda at harmonic n2
d2lam_n3=d2lam_n2_calc(d2lam_n1,n1,n3,phi);%based on guess, calc. d2lambda at harmonic n3
z_star_liq_n1=zqliq_n2_calc(z_star_liq_n,n,n1);%liq load impedance at harmonic n1
z_star_liq_n2=zqliq_n2_calc(z_star_liq_n,n,n2);%liq load impedance at harmonic n2
z_star_liq_n3=zqliq_n2_calc(z_star_liq_n,n,n3);%liq load impedance at harmonic n3
delf_n1_star_liq=delfnstar_liq_calc(f1,zq,z_star_liq_n1);%liq complex freq shift at harmonic n1
delf_n2_star_liq=delfnstar_liq_calc(f1,zq,z_star_liq_n2);%liq complex freq shift at harmonic n2
delf_n3_star_liq=delfnstar_liq_calc(f1,zq,z_star_liq_n3);%liq complex freq shift at harmonic n3
Rliq_n1=Rliq_calc(delf_n1_star_liq,n1,drho,f1,zq);%normalized complex freq shift at harmonic n1
Rliq_n2=Rliq_calc(delf_n2_star_liq,n2,drho,f1,zq);%normalized complex freq shift at harmonic n2
Rliq_n3=Rliq_calc(delf_n3_star_liq,n3,drho,f1,zq);%normalized complex freq shift at harmonic n3
Dn1=Dn_calc(d2lam_n1,phi);%Calculate the D coefficient at harmonic n1 (material property)
Dn2=Dn_calc(d2lam_n2,phi);%Calculate the D coefficient at harmonic n2 (material property)
Dn3=Dn_calc(d2lam_n3,phi);%Calculate the D coefficient at harmonic n3 (material property)
f_n1=f1*n1;%resonant frequency at n1 in Hz
f_n2=f1*n2;%resonant frequency at n2 in Hz
f_n3=f1*n3;%resonant frequency at n3 in Hz
%only used drho_n2 because drho_n1, n2, n3 should be similar
% harm_ratio_calc=real(master(d2lam_n2,phi,z_star_liq_n2,f_n2,drho))./real(master(d2lam_n1,phi,z_star_liq_n1,f_n1,drho));%calculate the harmonic ratio
% diss_ratio_calc=imag(master(d2lam_n3,phi,z_star_liq_n3,f_n3,drho))./real(master(d2lam_n3,phi,z_star_liq_n3,f_n3,drho));%calculate the dissipation ratio
harm_ratio_calc=real(master2(Dn2,Rliq_n2))./real(master2(Dn1,Rliq_n1));%calculate the harmonic ratio
diss_ratio_calc=imag(master2(Dn3,Rliq_n3))./real(master2(Dn3,Rliq_n3));%calculate the dissipation ratio
drho_calc1=drho_calc(handles.din.cursor.(['interp_harmfi',num2str(n1)]),zq,f1,n1,master2(Dn1,Rliq_n1));%calculate the drho value
F=[harm_ratio_calc-harm_ratio_exp;diss_ratio_calc-diss_ratio_exp;1000.*(drho_calc1-drho)];%calculates the difference between the calculated and experimental ratios

function handles=qcm_solver(handles,count)
%this function calculates the two variables, d2lambda and phi.
%handles: the gui handles structure
%count: keeps track which row to store the calculated viscoelasetic
%parameters in a variable
[radiotot,harmchoice]=radio_total(handles);%determine how many datasets and which harmonics that will be used to perform the calculation
f1=handles.din.constants.f1;%fundamentalr resonance freq in Hz
zq=handles.din.constants.zq;% the quartz load impedance
for dum0=1:length(radiotot)
    calc_table=nan(6,9);    %preallocate the calculation table
    handles.din.solved.harmchoice=harmchoice(dum0);
    [n1,n2,n3]=determine_harm_choice(harmchoice(dum0));%determine what harmonic datasets will be used to calculate viscoelastic parameters
    %calculate the liquid load impedance
    dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
    dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance at dgliq_harm
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at
    %extract out the guess parameters
    guess_table=get(handles.uitable4,'data');
    guess_values(1)=guess_table{4,2};%guess value for d2lambda at n1
    guess_values(2)=guess_table{2,2};%guess value for phi
    if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
        try
            drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
            drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
            drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
            guess_values(3)=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values and store it as aguess value for drho (kg/m^2)
            if isnan(guess_values(3))==0             
                guess_table{1,2}=guess_values(3)*1000;
                set(handles.uitable4,'data',guess_table);
            end%if isnan(drho_ave)==0
        catch
            guess_values(3)=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
        end%try
    else%if the manual guess for drho has been enabled
        guess_values(3)=guess_table{1,2}/1000;%use the drho value listed in the guess table instead
    end%if get(handles.edit_drho,'value')==0
       %calculate experimental harmonic and dissipation ratio
   try
        del_f_n2=handles.din.cursor.(['interp_harmfi',num2str(n2)]);%Df_n2
        del_f_n1=handles.din.cursor.(['interp_harmfi',num2str(n1)]);%Df_n1
        del_g_n3=handles.din.cursor.(['interp_harmgi',num2str(n3)]);%Dg_n3
        del_f_n3=handles.din.cursor.(['interp_harmfi',num2str(n3)]);%Df_n3
        harm_ratio_exp=(del_f_n2.*n1)./(del_f_n1.*n2);%experimental harmonic ratio
        diss_ratio_exp=del_g_n3./del_f_n3;%experimental dissipation ratio
        %calculate error in ratios
        [harm_ratio_error,diss_ratio_error,handles]=calc_ratio_error(handles,n1,n2,n3,del_f_n2,del_f_n1,del_g_n3,del_f_n3);
        %solve for phi and d2lambda
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp,diss_ratio_exp,x);
        options=handles.din.fit_options;
        [solved,resnorm,residual,exitflag,output,~,jacobian]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,options);
        output.resnorm=resnorm;% stores the squared 2-norm of the residue at x
        output.residual=residual;%stores the value of the residual at the solution x
        output.jacobian=jacobian;%store the jacobian matrix from fsolve
        output.exitflag=exitflag;%store the exitflag that was outputed from fsolve
        disp(['exitflag: ',num2str(exitflag),', resnorm: ',num2str(resnorm),', # of pts solved: ',num2str(count),', dataset: ',num2str(n1),num2str(n2),num2str(n3)]);
        disp(solved);
        if exitflag~=1&&resnorm>0.0001%stop the code if the exitflag is not 1 or if the sum of the squares of the residuals are greater than 1e-10
            set(handles.status,'string','Status: Exitflag is not 1!','backgroundcolor','y','foregroundcolor','r');
            disp('Exitflag is not 1!');
            return
        else
            set(handles.status,'string','Calculating...','backgroundcolor','k','foregroundcolor','r');
        end%if exitflag~=1
        handles.din.solved.output_misc_from_fsolve=output;%store the output structure in the handles structure
        handles.din.solved.solution=solved;%store the solution, solved(1): d2lam at n1 and solved(2): phi, solved(3): drho at n1
        %solve for the other viscoelastic parameters                
            unique_n=[n1;n2;n3;active_harm(handles)];
            unique_n=unique(unique_n,'stable');
        for dum=1:length(unique_n)            
            name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
            handles=calc_viscoelastic_par(handles,unique_n(dum),f1,z_star_liq_n,zq,solved,dgliq_harm,name,n1);%calculate the viscoelastic parameters
            
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
        %update the data tables int he gui
        set(handles.uitable2,'data',calc_table);%update the viscoelastic parameter table
        set(handles.uitable3,'data',pred_shifts);%update the predicted frequency shift table                        
    catch err_msg
        assignin('base','err_msg',err_msg);
        set(handles.status,'string','Status: Error in solving!','backgroundcolor','r','foregroundcolor','k');
        disp('Error in solving!');
        return
    end%try
    try%calculate the error
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp+harm_ratio_error,diss_ratio_exp+diss_ratio_error,x);
        [solved_e1,~,exitflage1,~,~,~]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp-harm_ratio_error,diss_ratio_exp-diss_ratio_error,x);
        [solved_e2,~,exitflage2,~,~,~]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp+harm_ratio_error,diss_ratio_exp-diss_ratio_error,x);
        [solved_e3,~,exitflage3,~,~,~]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp-harm_ratio_error,diss_ratio_exp+diss_ratio_error,x);
        [solved_e4,~,exitflage4,~,~,~]=lsqnonlin(fcns,guess_values,handles.din.lb,handles.din.ub,options);
        error0.solved_e1=solved_e1;
        error0.solved_e2=solved_e2;
        error0.solved_e3=solved_e3;
        error0.solved_e4=solved_e4;
        for dum0=1:4
            for dum3=1:length(unique_n)
                name=['n_',num2str(unique_n(dum3)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
                name_e=['solved_e',num2str(dum0)];
                handles=calc_viscoelastic_par(handles,unique_n(dum3),f1,z_star_liq_n,zq,error0.(name_e),dgliq_harm,name_e,n1);%calc the viscoelastic parameters
                if dum0~=1
                    %store the drho value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.drho(count,:)=[handles.din.cursor.selected_time,...
                        min([handles.din.stored_solutions.(name).error.drho(count,2),handles.din.solved.(name_e).drho]),...
                        max([handles.din.stored_solutions.(name).error.drho(count,3),handles.din.solved.(name_e).drho])];
                    %store the grho value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.grho(count,:)=[handles.din.cursor.selected_time,...
                        min([handles.din.stored_solutions.(name).error.grho(count,2),handles.din.solved.(name_e).grho]),...
                        max([handles.din.stored_solutions.(name).error.grho(count,3),handles.din.solved.(name_e).grho])];
                    %store the phi value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.phi(count,:)=[handles.din.cursor.selected_time,...
                        min([handles.din.stored_solutions.(name).error.phi(count,2),handles.din.solved.(name_e).phi]),...
                        max([handles.din.stored_solutions.(name).error.phi(count,3),handles.din.solved.(name_e).phi])];
                    %store the d2lamn value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.d2lam(count,:)=[handles.din.cursor.selected_time,...
                        min([handles.din.stored_solutions.(name).error.d2lam(count,2),handles.din.solved.(name_e).d2lam]),...
                        max([handles.din.stored_solutions.(name).error.d2lam(count,3),handles.din.solved.(name_e).d2lam])];
                    %store the complex pred. freq. values lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.complex_pred_freq(count,:)=[handles.din.cursor.selected_time,...
                        min([handles.din.stored_solutions.(name).error.complex_pred_freq(count,2),handles.din.solved.(name_e).complex_pred_freq]),...
                        max([handles.din.stored_solutions.(name).error.complex_pred_freq(count,3),handles.din.solved.(name_e).complex_pred_freq])];
                else
                    %store the drho lb and ub error and timepoint
                    handles.din.stored_solutions.(name).error.drho(count,:)=[handles.din.cursor.selected_time,...
                        handles.din.solved.(name_e).drho,handles.din.solved.(name_e).drho];
                    %store the grho lb and ub error and timepoint
                    handles.din.stored_solutions.(name).error.grho(count,:)=[handles.din.cursor.selected_time,...
                        handles.din.solved.(name_e).grho,handles.din.solved.(name_e).grho];
                    %store the phi value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.phi(count,:)=[handles.din.cursor.selected_time,...
                        handles.din.solved.(name_e).phi,handles.din.solved.(name_e).phi];
                    %store the d2lamn value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.d2lam(count,:)=[handles.din.cursor.selected_time,...
                        handles.din.solved.(name_e).d2lam,handles.din.solved.(name_e).d2lam];
                    %store the complex freq shift value lb and ub and timepoint
                    handles.din.stored_solutions.(name).error.complex_pred_freq(count,:)=[handles.din.cursor.selected_time,...
                        handles.din.solved.(name_e).complex_pred_freq,handles.din.solved.(name_e).complex_pred_freq];
                end%if count~=0
            end%for dum=1:length(unique_n)
        end%for dum0=1:length()
    catch err_msg
        assignin('base','err_msg',err_msg);
        set(handles.status,'string','Status: Error in calculating errors for harm and diss ratios!','backgroundcolor','r','foregroundcolor','k');
        disp('Error in calculating errors for harm and diss ratios!');
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

function handles=calc_viscoelastic_par(handles,n1,f1,z_star_liq_n,zq,solved,dgliq_harm,name,n0)
%this function calculates the viscoelastic parameters at harmonic n1 based on the solved
%variables, d2lam and phi
%handles: gui handles structure
%n1: the harmonic in which the parameters are calculated at
%f1: the fundamental resonance frequency
%z_star_liq_n: the calculated load impedance based on the dissipation shift
%between the bare xtal imersed in liquid and the bare xtal in air
%zq: load impedance of quartz
%solved: the solved parameters (d2lam and phi)
% drho_ave: the average (based on harmonics that were used to determine d2lam and phi) areal mass
%dgliq_harm: the harmonic in which the z_star_liq_n is associated
%name: the designated fieldname in which the calculations will be stored in
%the handles.din.solved structure
%n0: the harmonic associated with the solved d2lam
[nc1,nc2,nc3]=determine_harm_choice(handles.din.solved.harmchoice);%determine what harmonics are being used to calculate the viscoelastic parameters
handles.din.solved.(name).d2lam=d2lam_n2_calc(solved(1),n0,n1,solved(2));%store the d2lambda ratio at n1
handles.din.solved.(name).norm_delfstar=master(handles.din.solved.(name).d2lam,solved(2),zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1),f1*n1,solved(3));%normalized complex frequency shift at harmonic n1
if sum([nc1==n1,nc2==n1,nc3==n1])==0||handles.din.simulate==1%if n1 is not one of the harmonics that are being used to calculate the viscoeasltic parameters, use the calculated areal mass at nc1 to predict the frequency shift. 
    %This prevents the predicted freq shift from being the same as the exp. freq shift    
    name0=['n_',num2str(nc1),'_',num2str(nc1),num2str(nc2),num2str(nc3)];%nomenclature: n_<n1>_<harm used to calc properties>
    if handles.din.simulate==1&&isfield(handles.din.solved,name0)==1
        handles.din.solved.(name).drho=handles.din.solved.(name0).drho;%calculate drho for harmonic nc1 that was used in g/m^2    
    else        
        handles.din.solved.(name).drho=drho_calc(handles.din.cursor.(['interp_harmfi',num2str(n1)]),zq,f1,n1,handles.din.solved.(name).norm_delfstar).*1000;%calculate drho for each harmonic that was used in g/m^2    
    end
else 
    handles.din.solved.(name).drho=drho_calc(handles.din.cursor.(['interp_harmfi',num2str(n1)]),zq,f1,n1,handles.din.solved.(name).norm_delfstar).*1000;%calculate drho for each harmonic that was used in g/m^2    
end%if sum([nc1==n1,nc2==n1,nc3==n1])==0
handles.din.solved.(name).complex_pred_freq=handles.din.solved.(name).norm_delfstar.*sauerbrey(n1,f1,handles.din.solved.(name).drho./1000,zq);%calculate the complex frequency shifts (Hz)
handles.din.solved.(name).phi=solved(2);%store the viscoelastic phase angle at n1 (degrees)
handles.din.solved.(name).grho=grho_calc(handles.din.solved.(name).drho./1000,handles.din.solved.(name).d2lam,n1*f1,solved(2))./1000;%calculate grho at n1 in units of Pa-g/cm^3
handles.din.solved.(name).cf=zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)./(n1.*f1.*handles.din.solved.(name).drho);%calculate the "cf" term
handles.din.solved.(name).lam_rho=(sqrt(handles.din.solved.(name).grho.*1000)./(n1.*f1.*cosd(solved(2)/2))).*1000;%calculate \lambda\rho (um-g/cm^3)
handles.din.solved.(name).del2d=-1/imag(handles.din.solved.(name).d2lam.*2.*pi.*(1-1i.*tand(handles.din.solved.(name).phi./2)));


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
set(t,'columneditable',[true true],'celleditcallback',{@store_error,handles});

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
    try
        for dum=initial_time:solve_inc:extent
            handles=guidata(handles.figure1);
            timeline=handles.din.rawdata.freq_shift(:,1);
            timepoint0=timeline(dum)
            event.Position(1)=timepoint0;
            handles=rawfit(handles,harm,event);
            drawnow;
            h=guidata(findall(0,'type','figure','name','Refit panel'));
            h.guess_values_options.Value=4;
            handles.raw.index=find(harm_dataf(:,1)==timepoint0);
            guidata(handles.figure1,handles);            
            for repeat=1:5
                fit_button_callback(hObject,1,guidata(handles.figure1),harm,h.guess_values_options,h.a4,h.a1,h.a2,h.a3);
            end
            drawnow;
            handles=guidata(handles.figure1);
%             handles.raw.index=find(harm_dataf(:,1)==timepoint0);
            guidata(handles.figure1,handles);            
            accept_fcn(1,1,handles);
            drawnow;                                  
        end
        disp('Automated refitting process finished.');
    catch
    end
else
    disp('Raw spectras have not been loaded!');
end


% --------------------------------------------------------------------
function viewraw_ClickedCallback(hObject, eventdata, handles)
dcm=datacursormode(handles.figure1);
flag=0;
if strcmp(get(hObject,'state'),'on')
    handles.fit_all.Visible='on';
    try
        try        
            if exist([handles.raw.pathname,handles.raw.filename],'file')==2
                [~,filename,ext]=fileparts(handles.raw.filename);
                pathname=handles.raw.pathname;
            else
                err
            end
        catch        
            if exist(fileloc,'file')==2     
                [~,filename,ext]=fileparts(handles.raw.filename);            
                pathname=handles.din.filepath;
            else
                err
            end
        end
        disp('  ');
        disp('Registered raw file!');
    catch
        try
            while flag==0
                disp('Unable to locate raw file!');
                disp('Please select raw file to load');
                [filename,pathname,~]=uigetfile('.mat','Load raw data',handles.din.filepath);
                [~,filename,ext]=fileparts(filename);
                if strcmp(filename(end-12:end),'_raw_spectras')
                    flag=1;
                else
                    disp('Selected raw data file is not recognized!');
                    set(handles.status,'string','Status: Selected raw data file is not recognized',...
                        'backgroundcolor','r','foregroundcolor','k');
                end
            end
            set(handles.viewraw,'tooltipstring',['<html>view raw spectra <br/> filename: ',filename,'<br/> path: ',pathname,'</html>'],...
                'userdata',[{pathname},{filename}]);
            disp(['loading ',pathname,filename]);
            disp('Registered selected raw file!');
            set(handles.status,'string','Status: Registered selected raw file!','backgroundcolor','k','foregroundcolor','r');
        catch
            disp('Error loading data!');
            set(handles.status,'string','Status: Error loading data!','backgroundcolor','k','foregroundcolor','r');        
        end
    end    
    handles.raw.filename=[filename,'.mat'];
    handles.raw.pathname=pathname;
%     handles.spectra_file=matfile([pathname,filename]);
    guidata(handles.figure1,handles);
    set(dcm,'enable','on');
    set(dcm,'UpdateFcn',{@myupdatefcn,handles,hObject});
else
    set(dcm,'enable','off');
    set(handles.fit_all,'visible','off');
    delete(figure(1));
    delete(figure(2));
    delete(figure(3));
    r=findall(0,'type','hggroup');
    for dum=1:length(r)
        delete(r(dum));
    end
end

function txt=myupdatefcn(~,event,handles,viewraw)
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
if strcmp(get(viewraw,'state'),'on')
    handles=rawfit(handles,harm,event);
end
handles.raw.index=index;
guidata(handles.figure1,handles);

function handles=rawfit(handles,harm,event)
flag=0;
handles=guidata(handles.figure1);
h=guidata(figure(2));
delete(findall([figure(1) figure(2)],'userdata','mark'));
if isfield(h,'ver')==1    
    flag=1;
else
    clf(figure(1));set(figure(1),'units','normalized');pos=get(figure(1),'position');
    h.a1=axes;hold(h.a1,'on'); h.a2=axes;hold(h.a2,'on');
    xlabel(h.a1,'Frequency (Hz)','fontweight','bold','fontsize',14);
    ylabel(h.a1,'Conductance (mS)','fontweight','bold','fontsize',14,'color','b');
    xlabel(h.a2,'Frequency (Hz)','fontweight','bold','fontsize',14);
    dum=ylabel(h.a2,'Susceptance (mS)','fontweight','bold','fontsize',14,'color','r');
    set(dum,'rotation',-90,'verticalalignment','bottom');
end
[pathloc,fileloc,fileext]=fileparts(get(handles.filename_txt,'tooltipstring'));
try
    [~,filename,ext]=fileparts(handles.raw.filename);
    pathname=handles.raw.pathname;
    ind=find(handles.din.rawdata.freq_shift(:,1)==event.Position(1));
    timepoint=strrep(num2str(event.Position(1)),'.','dot');
    varname=sprintf([filename(1:end-13),'_t_%s_iq_1_ih_',num2str(0.5*(harm+1))],timepoint);
    disp(['Loading variable: ',varname]);
    tic
%     raw.(varname)=handles.spectra_file.(varname);
    raw=load([pathname,filename,ext],'-mat',varname);
    toc
    disp(['Variable loaded!']);
    set(handles.status,'string','Status: Variable loaded!','backgroundcolor','k','foregroundcolor','r');
    if flag==1
        set(findall(figure(1),'type','line','userdata','mark1'),'xdata',raw.(varname)(:,1),'ydata',raw.(varname)(:,2));
        set(findall(figure(1),'type','line','userdata','mark2'),'xdata',raw.(varname)(:,1),'ydata',raw.(varname)(:,3));
        set(findall(figure(2),'type','line','userdata','mark3'),'xdata',raw.(varname)(:,2),'ydata',raw.(varname)(:,3));
        set(figure(1),'name',[filename,'  View raw data for n=',num2str(harm),', ',num2str(event.Position(1)),' min',' index: ',num2str(ind)]);
    elseif flag==0
        plot(h.a1,raw.(varname)(:,1),raw.(varname)(:,2),'bx','userdata','mark1');
        plot(h.a2,raw.(varname)(:,1),raw.(varname)(:,3),'rx','userdata','mark2');  
        set(figure(1),'position',[pos(1) pos(2) 0.333 0.4],'name',...
            [filename,'  View raw data for n=',num2str(harm),', ',num2str(event.Position(1)),' min',' index: ',num2str(ind)],...
        'numbertitle','off');
        pos2=get(h.a1,'position');
        set(h.a1,'ycolor','b','position',[0.8*pos2(1) pos2(2) pos2(3) pos2(4)],'color','none');
        set(h.a2,'position',get(h.a1,'position'),'color','none','yaxislocation','right','ycolor','r');
        clf(figure(2));
        set(figure(2),'units','normalized','position',[pos(1)+0.333 pos(2) 0.165 0.4],'name','Refit panel','numbertitle','off');
        h.guess_values_options=uicontrol('parent',figure(2),'style','popupmenu',...
            'unit','normalized',...
            'fontweight','bold','fontsize',10,'string',[{'Gmax'};{'Derivative'};{'Bmax'};{'Previous values'};{'User-defined'}],...
            'horizontalalignment','left','tooltipstring','Choose guess values for curve fitting.',...
            'value',4,'position',[.025 0 0.45 0.08]);
        h.fit_button=uicontrol('parent',figure(2),'style','pushbutton','unit','normalized',...
            'position',[0.78 0.01 0.2 0.08],'string','Fit','fontweight','bold','tooltipstring',...
            'Fit the resonance peak');
        h.multi_button=uicontrol('parent',figure(2),'style','togglebutton','unit','normalized',...
            'position',[0.5 0.01 0.25 0.08],'string','Multi-peak','fontweight','bold','backgroundcolor',...
            [0.8 0.8 0.8]);           
        h.a3=axes;set(h.a3,'parent',figure(2),'position',[0.2 0.63 0.7 0.35],'fontsize',8);
        h.a4=axes;set(h.a4,'parent',figure(2),'position',[0.2 0.2 0.7 0.3],'fontsize',8);    
        xlabel(h.a3,'Conductance (mS)','fontweight','bold');
        ylabel(h.a3,'Susceptance (mS)','fontweight','bold');
        set(h.fit_button,'callback',{@fit_button_callback,handles,harm,h.guess_values_options,h.a4,h.a1,h.a2,h.a3});
        set(h.multi_button,'callback',{@multi,handles});
        plot(h.a3,raw.(varname)(:,2),raw.(varname)(:,3),'x','color',[0 0.5 0],'userdata','mark3');
    end        
    drawnow;    
    handles.raw.frequency=raw.(varname)(:,1);
    handles.raw.conductance=raw.(varname)(:,2);
    handles.raw.susceptance=raw.(varname)(:,3);
    h.ver='2layer_fit';
    guidata(findall(0,'type','figure','name','Refit panel'),h);
    guidata(handles.figure1,handles);
catch
    disp('Error in loading variable!');
    set(handles.status,'string','Status: Error in loading variable!','foregroundcolor','r','backgroundcolor','k');
end

function multi(hObject,~,handles)
if get(hObject,'value')==1
    handles=guidata(handles.figure1);
    clf(figure(3));set(figure(2),'units','normalized');pos0=get(figure(2),'position');
    set(figure(3),'toolbar','none','units','normalized','position',[pos0(1)-pos0(3) pos0(2)-0.23 pos0(3)*2 0.2],...
        'name','Multi-peak options','numbertitle','off','menubar','none');
    rnames={'1st','3rd','5th','7th','9th','11th'};
    cnames={'peak prominence sensitivity factor','peak minimum threshold','Max # peaks'};
    data=[handles.prefs.sensitivity(1:2:11);handles.prefs.peak_min(1:2:11);handles.prefs.num_peaks(1:2:11)];
    find_peak_options=uitable('units','normalized','position',[0.05 0.15 .9 .8],...
        'columnname',cnames,'rowname',rnames,'data',data','fontsize',10,...
        'columneditable',logical([1 1 1]),'celleditcallback',{@fpo,handles});%create the table
    uicontrol('style','pushbutton','string','table properties','units','normalized',...
        'position',[0.7 0.02 .25 .1],'callback',{@ins});
    set(figure(3),'closerequestfcn',{@fp_close,handles,hObject});
else
    delete(figure(3));
end

function fp_close(~,~,handles,multih)
set(multih,'value',0);
multi(multih,1,handles);

function fpo(hObject,~,handles)
data=get(hObject,'data');
handles.prefs.sensitivity(1:2:11)=data(:,1);
handles.prefs.peak_min(1:2:11)=data(:,2);
handles.prefs.num_peaks(1:2:11)=data(:,3);
guidata(handles.figure1,handles);
handles.prefs

function fit_button_callback(hObject,~,handles,harm,guess_values_options,a4,a1,a2,a3)
disp('Re-fitting the raw spectra!');
handles=guidata(handles.figure1);
handles.din.harmonic=harm;
freq=handles.raw.frequency;
conductance=handles.raw.conductance;
susceptance=handles.raw.susceptance;
combine_spectra=[freq,conductance,susceptance];
par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st'};
cla(a4);
del=findall(0,'userdata','mark');
for dum=1:length(del)
    delete(del(dum));
end
[G_fit,B_fit,G_l_sq,B_l_sq,combine_spectra,G_parameters,B_parameters,handles,I]...
    =Lorentzian_dynamic_fit(handles,freq,conductance,susceptance,combine_spectra,guess_values_options,a4);
if isempty(G_fit)==0    
    plot(a1,combine_spectra(:,1),G_fit,'k-','linewidth',2,'userdata','mark'); hold(a1,'on');    
    primary=lfun4_both_1([G_parameters(1:5),B_parameters(5)],combine_spectra(:,1));    
    offset=G_parameters(5);
    switch length(G_parameters)
        case 10
            offset=G_parameters(5)+G_parameters(10);
            plot(a1,G_parameters(6),G_parameters(9)+offset,'go','linewidth',1.5,'userdata','mark',...
                'markersize',12);
            second=lfun4_both_1([G_parameters(6:10),B_parameters(6:10)],combine_spectra(:,1));
            plot(a1,combine_spectra(:,1),second(:,1)+offset-G_parameters(10),'g--','linewidth',1.5,'userdata','mark');   
            primary=primary+G_parameters(10);
            pl2=text(G_parameters(6),G_parameters(9)+offset,{'2nd',' '});
            set(pl2,'edgecolor','none','backgroundcolor','none','horizontalalignment','center',...
            'verticalalignment','bottom','parent',a1,'color','g','fontweight','bold','userdata','mark');
            par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st',...
                'f0_2nd','gamma0_2nd','phi_2nd','Gmax_2nd','offset_2nd'};
        case 15
            offset=G_parameters(5)+G_parameters(10)+G_parameters(15);
            plot(a1,G_parameters(6),G_parameters(9)+offset,'go','linewidth',1.5,'userdata','mark',...
                'markersize',12);
            plot(a1,G_parameters(11),G_parameters(14)+offset,'go','linewidth',1.5,'userdata','mark',...
                'markersize',12);
            second=lfun4_both_1([G_parameters(6:10),B_parameters(6:10)],combine_spectra(:,1));
            plot(a1,combine_spectra(:,1),second(:,1)+offset-G_parameters(10),'g--','linewidth',1.5,'userdata','mark');
            pl2=text(G_parameters(6),G_parameters(9)+offset,{'2nd',' '});
            set(pl2,'edgecolor','none','backgroundcolor','none','horizontalalignment','center',...
            'verticalalignment','bottom','parent',a1,'color','g','fontweight','bold','userdata','mark');
            third=lfun4_both_1([G_parameters(11:15),B_parameters(11:15)],combine_spectra(:,1));
            plot(a1,combine_spectra(:,1),third(:,1)+offset-G_parameters(15),'g--','linewidth',1.5,'userdata','mark');
            pl3=text(G_parameters(11),G_parameters(14)+offset,{'3rd',' '});
            set(pl3,'edgecolor','none','backgroundcolor','none','horizontalalignment','center',...
            'verticalalignment','bottom','parent',a1,'color','g','fontweight','bold','userdata','mark');
            primary=primary+G_parameters(10)+G_parameters(15);
            par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st',...
            'f0_2nd','gamma0_2nd','phi_2nd','Gmax_2nd','offset_2nd',...
            'f0_3rd','gamma0_3rd','phi_3rd','Gmax_3rd','offset_3rd'};
    end%switch
    plot(a1,G_parameters(1),G_parameters(4)+offset,'mx','linewidth',1.5,'userdata','mark',...
        'markersize',12);
    pl1=text(G_parameters(1),G_parameters(4)+offset,{'1st',' '},'userdata','mark');
    set(pl1,'edgecolor','none','backgroundcolor','none','horizontalalignment','center',...
        'verticalalignment','bottom','parent',a1,'color','m','fontweight','bold','userdata','mark');
    plot(a1,combine_spectra(:,1),primary(:,1),'m--','linewidth',1.5,'userdata','mark');
    plot(a2,combine_spectra(:,1),B_fit,'k-','linewidth',2,'userdata','mark');
    hold(a3,'on');
    plot(a3,G_fit,B_fit,'k-','userdata','mark');
    linkaxes([a1 a2],'x');
    uistack(a1,'up');
    bare_ref=get(handles.bare_table,'data');
    ref=[bare_ref(0.5*(harm+1),1),bare_ref(0.5*(harm+1),2)];
    harm_dataf=handles.din.(['harmf',num2str(harm)]);
    harm_datag=handles.din.(['harmg',num2str(harm)]);    
    str={['\bf f: ',num2str(G_parameters(1)),' Hz'],...
        ['\bf\Deltaf: ',num2str(G_parameters(1)-ref(1)),' Hz'],...
        ['\bf\Deltaf_c_o_r_r_e_c_t_i_o_n: ',...
        num2str((G_parameters(1)-ref(1))-harm_dataf(handles.raw.index,2)),' Hz'],...
        ['\bf\Gamma: ',num2str(G_parameters(2)),' Hz'],...
        ['\bf\Delta\Gamma: ',num2str(G_parameters(2)-ref(2)),' Hz'],...
        ['\bf\Delta\Gamma_c_o_r_r_e_c_t_i_o_n: ',...
        num2str((G_parameters(2)-ref(2))-harm_datag(handles.raw.index,2)),' Hz']};
    c=get(a2,'children');d=findobj(c,'type','text');
    uistack(a1,'up');
    if isempty(d)==0
        delete(d);
    end
    txt=text(1,1,' ');
    set(txt,'parent',a2,'units','normalized','position',[0.7 0.9 0],'string',str,'userdata','mark');
    handles.raw.current_Gparameters=G_parameters;
    handles.raw.current_Bparameters=B_parameters;
    handles.raw.Glsq=G_l_sq;
    handles.raw.Blsq=B_l_sq;
    guidata(handles.figure1,handles);
    assignin('base','G_parameters',G_parameters);
    assignin('base','B_parameters',B_parameters);
    assignin('base','Parameter_labels',par_labels);
    disp('Fitting parameters have been exported to the workspace.');
    disp('Current fitting parameters (row 1: G_par, row2: B_par):');
    disp(par_labels);
    disp([G_parameters;B_parameters]);
    disp(' ');
    accept=uicontrol(get(a1,'parent'),'style','pushbutton','userdata','mark',...
        'string','Accept fit','callback',{@accept_fcn,handles});
    del=uicontrol(get(a1,'parent'),'style','pushbutton','userdata','mark',...
        'string','Delete datapoint','callback',{@del_fcn,handles});
    accept.Position(1)=5; accept.Position(2)=5;
    del.Position(1)=70; del.Position(2)=5; del.Position(3)=del.Position(3)*1.5;
    h.accept=accept;
    h.del=del;
    guidata(accept.Parent,h);
end

function accept_fcn(hObject,event,handles)
harm=handles.din.harmonic;%current harmonic
ref=handles.din.rawdata.freq_shift_ref(:,(harm+1)/2);
fdata=handles.din.(['harmf',num2str(harm)]);%extract current freq shift data
gdata=handles.din.(['harmg',num2str(harm)]);%extract current gamma shift data
index=handles.raw.index;%extract out index of datapoint that will be modified
ref_table=get(handles.bare_table,'data');%get ref. data
ref=ref_table(0.5*(harm+1),:);
Gpar=handles.raw.current_Gparameters;%new fitted G parameters
Bpar=handles.raw.current_Bparameters;%new fitted B parameters
Glsq=handles.raw.Glsq;
Blsq=handles.raw.Blsq;
fdata(index,2)=nanmean([Gpar(1),Bpar(1)])-ref(1);%update delta freq
gdata(index,2)=nanmean([Gpar(2),Bpar(2)])-ref(2);%update delta gamma
%update the handles structure
handles.din.(['harmf',num2str(harm)])=fdata;
handles.din.(['harmg',num2str(harm)])=gdata;
index2=find(handles.din.rawdata.freq_shift(:,1)==fdata(index(1)));
handles.din.rawdata.freq_shift(index2,harm+1)=fdata(index,2);
handles.din.rawdata.freq_shift(index2,harm+2)=gdata(index,2);
handles.din.rawdata.abs_freq(index2,harm+1)=fdata(index,2)+ref(1);
handles.din.rawdata.abs_freq(index2,harm+2)=gdata(index,2)+ref(2);
handles.din.rawdata.chisq_values(index2,harm+1)=sum(Glsq);
handles.din.rawdata.chisq_values(index2,harm+2)=sum(Blsq);
guidata(handles.figure1,handles);%update the handles structure
%save the changes into the imported .mat file
pathname=handles.din.filepath;
filename=handles.din.filename;
freq_shift=handles.din.rawdata.freq_shift;
abs_freq=handles.din.rawdata.abs_freq;
chisq_values=handles.din.rawdata.chisq_values;
freq_shift_ref=handles.din.rawdata.freq_shift_ref;
reference=handles.din.rawdata.reference;
version=handles.din.rawdata.version;
save([pathname,filename],'freq_shift','abs_freq','chisq_values','freq_shift_ref','reference','version');
% export=matfile([pathname,filename],'writable',true);
% export.freq_shift=handles.din.rawdata.freq_shift;
% export.abs_freq=handles.din.rawdata.abs_freq;
% export.chisq_values=handles.din.rawdata.chisq_values;
% export.(['mod_',datestr(now,'yyyy_mm_dd_HH_MM_SS')])={['Modifications made on datapoint: ',num2str(index)],...
%     ['Modification performed on: ',date],...
%     ['Program versions: ',get(handles.text4,'string')]};
disp('Dataset and file updated');
disp('Updating plots');
myharmfcn(handles.(['harm',num2str(harm)]),handles);
disp('Plots updated')

function del_fcn(hObject,evemt,handles)
harm=handles.din.harmonic;%current harmonic
fdata=handles.din.(['harmf',num2str(harm)]);%extract current freq shift data
gdata=handles.din.(['harmg',num2str(harm)]);%extract current gamma shift data
index=handles.raw.index;%extract out index of datapoint that will be modified
fdata(index,2)=nan;%update delta freq
gdata(index,2)=nan;%update delta gamma
%update the handles structure
handles.din.(['harmf',num2str(harm)])=nan;
handles.din.(['harmg',num2str(harm)])=nan;
index2=find(handles.din.rawdata.freq_shift(:,1)==fdata(index(1)));
handles.din.rawdata.freq_shift(index2,harm+1)=fdata(index,2);
handles.din.rawdata.freq_shift(index2,harm+2)=gdata(index,2);
handles.din.rawdata.abs_freq(index2,harm+1)=nan;
handles.din.rawdata.abs_freq(index2,harm+2)=nan;
handles.din.rawdata.chisq_values(index2,harm+1)=nan;
handles.din.rawdata.chisq_values(index2,harm+2)=nan;
guidata(handles.figure1,handles);%update the handles structure
pathname=handles.din.filepath;
filename=handles.din.filename;
%save the changes into the imported .mat file
pathname=handles.din.filepath;
filename=handles.din.filename;
export=matfile([pathname,filename],'writable',true);
export.freq_shift=handles.din.rawdata.freq_shift;
export.abs_freq=handles.din.rawdata.abs_freq;
export.chisq_values=handles.din.rawdata.chisq_values;
export.refit_log={['Modifications made on datapoint: ',num2str(index)],...
    ['Modification performed on: ',date],...
    ['Program versions: ',get(handles.text4,'string')]};
disp('Dataset and file updated');
disp('Updating plots');
myharmfcn(handles.(['harm',num2str(harm)]),handles);
disp('Plots updated')
delete(figure(2));
delete(figure(1));

function [G_fit,B_fit,G_l_sq,B_l_sq,combine_spectra,G_parameters,B_parameters,handles,I]=Lorentzian_dynamic_fit(handles,freq,conductance,susceptance,combine_spectra,guess_values_options,a4)
factor_range_fit=handles.din.fit_factor_range;
%make sure to change user-defined values to previous values if the user is
%conducting the measurements. This prevents the program from pausing and
%waiting for the user to input the guess parameters in the middle of a
%measurement.
G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;  I=1;
switch get(guess_values_options,'value')
    case 1%Guess value based on max conductance
        [guess,f0,gamma0]=G_guess(freq,conductance,susceptance,handles,'Conductance (mS)',a4);
        if isempty(guess)==1
            G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
            return
        end%if isempty(guess)
        I=find(freq>=(f0-gamma0*factor_range_fit)&freq<=(gamma0*factor_range_fit+f0)); 
        try %fitting with Gmax initial guesses
            [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
            combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
        catch% tryGuess values based on the Derivative of the Fit
            disp('Fitting based on the Gmax guess  failed!!');
            disp('Attempting to use derivative values to fit...');
            [p,freq_mod,modulus,~,~]=deriv_guess(freq,conductance,susceptance,handles,a4);
            if isempty(p)==1
                G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
                return
            end%if isempty(guess)
            I=find(freq_mod>=(f0-gamma0*factor_range_fit)&freq_mod<=(gamma0*factor_range_fit+f0)); 
            [~,~,test]=fit_spectra_con(p,freq_mod,modulus,I,1);
            guess=[test(1) test(2) p(3:4) mean([conductance(1) conductance(end)])];%guess values
            try% tryGuess values based on the Derivative of the Fit
                [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
                combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                disp('Gmax guess values suceeded!');
            catch%if fit fails, output nan arrays
                disp('Fit failed!');
                G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;                                                
            end%try
        end%try
    case 2  %Guess values based on the Derivative of the Fit
        [guess,freq_mod,modulus,f0,gamma0]=deriv_guess(freq,conductance,susceptance,handles,a4);
        if isempty(guess)==1
            G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
            return
        end%if isempty(guess)
        I=find(freq_mod>=(f0-gamma0*factor_range_fit)&freq_mod<=(gamma0*factor_range_fit+f0)); 
        try
            [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq conductance susceptance],guess,I);
            combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
        catch%try GMAx as guess values
            disp('Fitting based on the derivative failed!!');
            disp('Attempting to use Gmax guess values to fit...');
            [guess,~,~]=G_guess(freq,conductance,susceptance,handles,'Conductance (mS)',a4);       
            if isempty(guess)==1
                G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
                return
            end%if isempty(guess)
            try
                [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq conductance susceptance],guess,I);
                combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                disp('Gmax guess values suceeded!');
            catch%if fit fails, output nan arrays
                disp('Fit failed!');
                G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;
            end%try            
        end%try           
    case 4%Guess value base on the previous fit values
        if isfield(handles.raw,'current_Gparameters')==1&&isfield(handles.raw,'current_Bparameters')
            disp('Previous guess parameters found.');
            G_prev=handles.raw.current_Gparameters;
            prev_par=[handles.raw.current_Gparameters;handles.raw.current_Bparameters];
            switch length(prev_par)
                case 5
                    par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st'};
                case 10
                    par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st',...
                    'f0_2nd','gamma0_2nd','phi_2nd','Gmax_2nd','offset_2nd'};
                case 15
                    par_labels={'f0_1st','gamma0_1st','phi_1st','Gmax_1st','offset_1st',...
                    'f0_2nd','gamma0_2nd','phi_2nd','Gmax_2nd','offset_2nd',...
                    'f0_3rd','gamma0_3rd','phi_3rd','Gmax_3rd','offset_3rd'};
            end
            disp('Previous fitting parameters (row 1: G_par, row2: B_par):');
            disp(par_labels)
            disp(prev_par);
            guess=mean(prev_par,1);
            if handles.prefs.simul_peak==1
                if length(guess)==5
                    guess=[guess prev_par(2,5)];
                elseif length(guess)==10
                    guess=[guess prev_par(2,5) prev_par(2,10)];
                elseif length(guess)==15
                    guess=[guess prev_par(2,5) prev_par(2,10) prev_par(2,15)];
                end% if length(guess)==5
            end%if handles.prefs.simul_peak==1
            I=find(freq>=(G_prev(1)-G_prev(2)*factor_range_fit)&freq<=(G_prev(2)*factor_range_fit+G_prev(1))); 
            try
                [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq conductance susceptance],guess,I);
            catch
                [guess,freq_mod,modulus,f0,gamma0]=deriv_guess(freq,conductance,susceptance,handles,a4);
                if isempty(guess)==1
                    return
                end%if isempty(guess)
                I=find(freq_mod>=(f0-gamma0*factor_range_fit)&freq_mod<=(gamma0*factor_range_fit+f0)); 
                try
                    [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq conductance susceptance],guess,I);
                    combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                catch%if fit fails, output nan arrays
                    disp('Fit failed!');
                    G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                    B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;
                end%try  
            end
            combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
        else%if not try to fit the data by choosing guess values from the derivative of the polar plot
            disp('Previous guess values not found. Guess values will be chosen from the derivative of the polar plot');
            [guess,freq_mod,modulus,f0,gamma0]=deriv_guess(freq,conductance,susceptance,handles,a4);
            if isempty(guess)==1
                return
            end%if isempty(guess)
            I=find(freq_mod>=(f0-gamma0*factor_range_fit)&freq_mod<=(gamma0*factor_range_fit+f0)); 
            try
                [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq conductance susceptance],guess,I);
                combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
            catch%if fit fails, output nan arrays
                disp('Fit failed!');
                G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;
            end%try             
        end%if isfield(handles.din,'G_prev')
        if isempty(guess)==1
            G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
            return
        end%if isempty(guess)
    case 3%Use the susceptance spectra to find the guess values
        [guess,f0,gamma0]=G_guess(freq,susceptance,conductance,handles,'Susceptance (mS)',a4);
        if isempty(guess)==1
            G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
            return
        end%if isempty(guess)
        I=find(freq>=(f0-gamma0*factor_range_fit)&freq<=(gamma0*factor_range_fit+f0)); 
        try %fitting with Gmax initial guesses
            [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
            combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
        catch% tryGuess values based on the Derivative of the Fit
            disp('Fitting based on the Bmax guess  failed!!');
            disp('Attempting to use derivative values to fit...');
            [p,freq_mod,modulus,~,~]=deriv_guess(freq,conductance,susceptance,handles,a4);
            if isempty(p)==1
                G_fit=[];B_fit=[];G_l_sq=[];B_l_sq=[];combine_spectra=[];G_parameters=[];B_parameters=[];I=[];
                return
            end%if isempty(guess)
            [~,~,test]=fit_spectra_con(p,freq_mod,modulus,1);
            guess=[test(1) test(2) p(3:4) mean([conductance(1) conductance(end)])];%guess values
            try% tryGuess values based on the Derivative of the Fit
                [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
                combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                disp('Gmax guess values suceeded!');
            catch%if fit fails, output nan arrays
                disp('Fit failed!');
                G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;                                                
            end%try
        end%try
    case 5%Guess values based on user defined parameters
        user_peaks=questdlg('How many peaks do you wish to define?','User defined peaks',...
            'Two','Three','Cancel','Two');%ask how many peaks to fit
        switch user_peaks
            case 'Two'%fit 2 peaks
                %provide starting guess values for user to edit
                try
                    [peak_detect,index]=findpeaks(smooth(smooth(conductance),7),'minpeakprominence',peak_sens((handles.din.harmonic+1)/2,1),...
                        'minpeakheight',(max(conductance)-min(conductance))*peak_sens((handles.din.harmonic+1)/2,2)+min(conductance));
                catch
                    [peak_detect,index]=findpeaks(smooth(conductance),'sortstr','descend');
                end%try 
                Gmax=peak_detect(1);%find peak of curve, Gmax refers to the maximum of the conductance
                f0=freq(index(1));%finds freq at which Gmax happens
                halfg=(Gmax-min(conductance))./2+min(conductance);%half of the Gmax
                halfg_freq=freq(find(abs(halfg-conductance)==min(abs((halfg-conductance))),1));%estimate of Gamma value
                gamma0=abs(halfg_freq-f0);%Guess for gamma, HMHW of peak
                offset=0; phi=0;
                try
                    guess=get(handles.(['X',num2str(handles.din.harmonic)]),'userdata');
                    guess=[[guess(1,1:5)';guess(2,5)],[guess(1,6:10)';guess(2,10)]];
                catch                                       
                    try
                        guess=[[f0;gamma0;phi;Gmax;offset;offset],...
                        [freq(index(2));gamma0/4;phi;peak_detect(2);offset;offset]];           
                    catch
                        guess=[[f0;gamma0;phi;Gmax;offset;offset],...
                        [freq(index(1));gamma0/4;phi;peak_detect(1);offset;offset]];     
                    end%try
                end
                figure(997);clf(figure(997));pos=get(figure(997),'position');set(figure(997),'position',[pos(1) pos(2) 400 300]);
                ud_table=uitable(figure(997),'units','normalized','position',[0.05 0.05 0.9 0.9],'data',guess,...
                    'RowName',[{'freq.'},{'gamma0'},{'\phi'},{'Gmax'},{'Cond. offset'},{'Sus. Offset'}],...
                    'ColumnName',[{'Peak 1'},{'Peak 2'}],'CellEditCallback',{@ud_values,handles},...
                    'ColumnEditable',[true true]);%create table for the user to input the guess values
                ud_confirm=uicontrol('style','pushbutton','string','OK','parent',figure(997),'backgroundcolor',[0.75 0.75 0.75],...
                    'callback',@ud_confirm_callback,'userdata',0);
                waitfor(ud_confirm,'userdata');%wait for the user until the "OK" button is pressed
                set(figure(997),'visible','off');drawnow;
                guess=get(ud_table,'data');%extract the user-defined guess values
                disp('Using user-defined values');
                I=find(freq>=(f0-gamma0*factor_range_fit)&freq<=(gamma0*factor_range_fit+f0)); 
                guess=[guess(1,1),guess(2,1),guess(3,1),guess(4,1),guess(5,1),guess(1,2),...
                    guess(2,2),guess(3,2),guess(4,2),guess(5,2),guess(6,1),guess(6,2)];%reformat the guess array
                %set the number of peaks to fit to 2
                handles.prefs.num_peaks(handles.din.harmonic)=2;
                guidata(handles.figure1,handles);
                try% tryGuess values based on user-input values
                    [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
                    combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                    disp('User-defined fit completed!');
                catch%if fit fails, output nan arrays
                    disp('Fit failed!');
                    G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                    B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;                                    
                end%try
            case 'Three'%fit 3 peaks
                %provide starting guess values for user to edit
                try
                    [peak_detect,index]=findpeaks(smooth(smooth(conductance),7),'minpeakprominence',peak_sens((handles.din.harmonic+1)/2,1),...
                        'minpeakheight',(max(conductance)-min(conductance))*peak_sens((handles.din.harmonic+1)/2,2)+min(conductance));
                catch
                    [peak_detect,index]=findpeaks(smooth(conductance),'sortstr','descend');
                end%try 
                Gmax=peak_detect(1);%find peak of curve
                f0=freq(index(1));%finds freq at which Gmax happens
                halfg=(Gmax-min(conductance))./2+min(conductance);%half of the Gmax
                halfg_freq=freq(find(abs(halfg-conductance)==min(abs((halfg-conductance))),1));
                gamma0=abs(halfg_freq-f0);%Guess for gamma, HMHW of peak
                offset=0; phi=0;
                try
                    guess=get(handles.(['X',num2str(handles.din.harmonic)]),'userdata');
                    guess=[[guess(1,1:5)';guess(2,5)],[guess(1,6:10)';guess(2,10)],[guess(1,11:15)';guess(2,15)]];
                catch                                       
                    try
                        guess=[[f0;gamma0;phi;Gmax;offset;offset],...
                        [freq(index(2));gamma0/4;phi;peak_detect(2);offset;offset],...
                        [freq(index(3));gamma0/4;phi;peak_detect(3);offset;offset]];           
                    catch
                        guess=[[f0;gamma0;phi;Gmax;offset;offset],...
                        [freq(index(1));gamma0/4;phi;peak_detect(1);offset;offset],...
                        [freq(index(1));gamma0/4;phi;peak_detect(1);offset;offset]];     
                    end%try
                end       
                figure(997);clf(figure(997));pos=get(figure(997),'position');set(figure(997),'position',[pos(1) pos(2) 450 300]);
                ud_table=uitable(figure(997),'units','normalized','position',[0.05 0.05 0.9 0.9],'data',guess,...
                    'RowName',[{'freq.'},{'gamma0'},{'\phi'},{'Gmax'},{'Cond. offset'},{'Sus. Offset'}],...
                    'ColumnName',[{'Peak 1'},{'Peak 2'},{'Peak 3'}],'CellEditCallback',{@ud_values,handles},...
                    'ColumnEditable',[true true true]);%create table for the user to input the guess values
                ud_confirm=uicontrol('style','pushbutton','string','OK','parent',figure(997),'backgroundcolor',[0.75 0.75 0.75],...
                    'callback',@ud_confirm_callback,'userdata',0);
                waitfor(ud_confirm,'userdata');%wait for the user until the "OK" button is pressed
                set(figure(997),'visible','off');drawnow;
                guess=get(ud_table,'data');%extract the user-defined guess values
                disp('Using user-defined values');
                I=find(freq>=(f0-gamma0*factor_range_fit)&freq<=(gamma0*factor_range_fit+f0)); 
                guess=[guess(1,1),guess(2,1),guess(3,1),guess(4,1),guess(5,1),guess(1,2),...
                    guess(2,2),guess(3,2),guess(4,2),guess(5,2),...
                    guess(1,3),guess(2,3),guess(3,3),guess(4,3),guess(5,3),guess(6,1),guess(6,2),guess(6,3)];%reformat the guess array
                %set the number of peaks to fit to 3                
                handles.prefs.num_peaks(handles.din.harmonic)=3;
                guidata(handles.figure1,handles);
                try% tryGuess values based on user-input values
                    [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,[freq,conductance,susceptance],guess,I);
                    combine_spectra=[freq,conductance,susceptance,G_fit,B_fit,G_l_sq,B_l_sq];%put everything in one variable
                    disp('User-defined fit completed!');
                catch%if fit fails, output nan arrays
                    disp('Fit failed!');
                    G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
                    B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;                                    
                end%try
            case 'Cancel'%cancel user defined values and bring back to default Gmax choice
                G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;%output nan arrays
                B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;  I=1;
        end%switch user_peaks
end%switch 
if sum(isnan(G_parameters))==0&&sum(isnan(B_parameters))>0%This if statemens ensures that there will not be dimenaional mismatch
    B_parameters=nan(1,length(G_parameters));
end%if sum(isnan(G_parameters))==0&&sum(isnan(B_parameters))>0


function ud_values(hObject,event,handles)

function ud_confirm_callback(hObject,~)
set(hObject,'userdata',randi(10));

function [guess,f0,gamma0]=G_guess(freq,conductance,susceptance,handles,ylab,a4)
%this function finds the guess values based on the max conductance value of
%the raw spectra
phi=0;%Assume rotation angle is 0
offset=0;%Assume offset value is 0
peak_sens=handles.prefs.num_peaks;
peak_sens2=handles.prefs.sensitivity;
peak_sens3=handles.prefs.peak_min;
try
    [peak_detect,index]=findpeaks(smooth(smooth(conductance),7),'minpeakprominence',peak_sens2(handles.din.harmonic),...
        'minpeakheight',(max(conductance)-min(conductance))*peak_sens3(handles.din.harmonic)+min(conductance));
catch
    [peak_detect,index]=findpeaks(smooth(conductance),'sortstr','descend');
end%try
flag=preview_peak_identification(freq,conductance,index,ylab,handles,a4);
if flag==1
    guess=[]; f0=[]; gamma0=[];
    return
end%flag==1
Gmax=peak_detect(1);%find peak of curve
f0=freq(index(1));%finds freq at which Gmax happens
halfg=(Gmax-min(conductance))./2+min(conductance);%half of the Gmax
halfg_freq=freq(find(abs(halfg-conductance)==min(abs((halfg-conductance))),1));
gamma0=abs(halfg_freq-f0);%Guess for gamma, HMHW of peak
guess=[f0 gamma0 phi Gmax offset];%consolidate the guess values into a single variable
if peak_sens(handles.din.harmonic)>=1&&size(peak_detect,1)>0
    guess=[guess 0];
end%if peak_sens((handles.din.harmonic+1)/2,3)==1&&size(peak_detect,1)==1
if peak_sens(handles.din.harmonic)>=2&&size(peak_detect,1)>1
    guess=[f0 gamma0 phi Gmax offset...
        freq(index(2)) gamma0/4 phi peak_detect(2) offset offset offset];
end%if peak_sens((handles.din.harmonic+1)/2,3)>=2&&size(peak_detect,1)>1
if peak_sens(handles.din.harmonic)>=3&&size(peak_detect,1)>2
    guess=[f0 gamma0 phi Gmax offset...
        freq(index(2)) gamma0/2 phi peak_detect(2) offset...
        freq(index(3)) gamma0/2 phi peak_detect(3) offset offset offset offset];
end%if peak_sens((handles.din.harmonic+1)/2,3)>=3&&size(peak_detect,1)>2

function [guess,freq_mod,modulus,f0,gamma0]=deriv_guess(freq,conductance,susceptance,handles,a4)
%this function finds the guess values based on the derivative of the raw
%spectra
phi=0;%Assume rotation angle is 0
offset=0;%Assume offset value is 0
modulus=sqrt((diff(conductance)).^2+(diff(susceptance)).^2);
freq_mod=freq(1:end-1)+diff(freq)./2;
peak_sens=handles.prefs.num_peaks;
peak_sens2=handles.prefs.sensitivity;
peak_sens3=handles.prefs.peak_min;
try
    [peak_detect,index]=findpeaks(smooth(modulus,9),'minpeakprominence',peak_sens2(handles.din.harmonic),...
        'minpeakheight',(max(smooth(modulus,9))-min(smooth(modulus,9)))*peak_sens3(handles.din.harmonic)+min(smooth(modulus,9)));
catch
    [peak_detect,index]=findpeaks(smooth(modulus,9),'sortstr','descend');
end%try
flag=preview_peak_identification(freq_mod,modulus,index,'Modulus (mS)',handles,a4);
if flag==1
    guess=[]; f0=[]; gamma0=[];freq_mod=[];modulus=[];
    return
end%flag==1
if isempty(peak_detect)
    disp('No peak detected!');
    guess=[]; f0=[]; gamma0=[];freq_mod=[];modulus=[];
    set(handles.status,'string','Status: No peak detected!','foregroundcolor','k','backgroundcolor','r');
    return
end%if isempty(peak_detect)
modulus=smooth(modulus,9);%smooth out the dataset
modulus_max=peak_detect(1);%find peak of curve
f0=freq_mod(index(1));%finds freq at which Gmax happens
halfg=(modulus_max-min(modulus))./2+min(modulus);%half of the Gmax
halfg_freq=freq_mod(find(abs(halfg-modulus)==min(abs((halfg-modulus))),1));
gamma0=abs(halfg_freq-f0);%Guess for gamma, HMHW of peak
phi=asind(conductance(1)/(sqrt((conductance(1))^2+(susceptance(1))^2)));%guess of the phase angle between the conductance and susceptance
guess=[f0 gamma0 phi modulus_max offset];%consolidate the guess values into a single variable
if peak_sens(handles.din.harmonic)>=1&&size(peak_detect,1)>0
    guess=[guess 0];
end%if peak_sens((handles.din.harmonic+1)/2,3)==1&&size(peak_detect,1)==1
if peak_sens(handles.din.harmonic)>=2&&size(peak_detect,1)>1
    guess=[f0 gamma0 phi modulus_max offset...
        freq(index(2)) gamma0/2 phi peak_detect(2) offset offset offset];
end%if peak_sens((handles.din.harmonic+1)/2,3)==2&&size(peak_detect,1)>1
if peak_sens(handles.din.harmonic)>=3&&size(peak_detect,1)>2
    guess=[f0 gamma0 phi modulus_max offset...
        freq(index(2)) gamma0/2 phi peak_detect(2) offset...
        freq(index(3)) gamma0/2 phi peak_detect(3) offset offset offset offset];
end%if peak_sens((handles.din.harmonic+1)/2,3)>=3&&size(peak_detect,1)>2

function flag=preview_peak_identification(freq,ydata,index,ylabel_str,handles,a4)
flag=0;
plot(a4,freq,smooth(smooth(ydata,7)));hold on;
plot(a4,freq(index),ydata(index),'-x','linewidth',2,'markersize',10);
xlabel('Freq (Hz)','fontweight','bold');
ylabel(ylabel_str,'fontweight','bold');
L=legend('Smoothed data',['Identified peaks: ',num2str(length(index))],'location','best');
set(L,'fontsize',8,'color','none')
drawnow;
button=questdlg('Continue or cancel fitting?','Fitting paused...','Continue','Cancel','Continue');
switch button
    case 'Continue'
        flag=0;
    case 'Cancel'
        flag=1;
        disp('Peak fitting canceled!');
end

function [G_fit,G_parameters,G_l_sq,B_fit,B_parameters,B_l_sq]=fit_spectra(handles,raw_data,guess,I)
%This function tries to fit the raw spectra using the provided guess values
freq=raw_data(:,1);
conductance=raw_data(:,2);
susceptance=raw_data(:,3);
if handles.prefs.simul_peak==0
    [G_fit,G_residual,G_parameters]=fit_spectra_con(guess,freq,conductance,I,1);            
    G_l_sq=(G_residual.^2)./((str2double(get(handles.num_datapoints,'string'))-1));%chi-squared calculation
    if get(handles.fit_B_radio,'value')==1
        [B_fit,B_residual,B_parameters]=fit_spectra_sus(guess,freq,susceptance,I,handles.prefs.show_GB);
        B_l_sq=(B_residual.^2)./((str2double(get(handles.num_datapoints,'string'))-1));%chi-squared calculation
    else%set susceptance variables to nan values
        B_fit=NaN(size(G_fit,1),size(G_fit,2));
        B_l_sq=B_fit;
        B_parameters=[NaN NaN NaN NaN NaN];
    end%if get(handles.fit_B_radio,'value')==1
elseif handles.prefs.simul_peak==1%run this block of code if the fitting G and B simultaneously option is turned on
    tic
    [GB_fit,GB_residual,GB_parameters]=fit_spectra_both(guess,freq,conductance,susceptance,handles.prefs.num_peaks,I,handles);
    toc
    try
        if length(GB_parameters)==6
            G_parameters=GB_parameters(1:5);        B_parameters=[GB_parameters(1:4) GB_parameters(6)];
        elseif length(GB_parameters)==12
            G_parameters=GB_parameters(1:10);       B_parameters=[GB_parameters(1:4) GB_parameters(11) GB_parameters(6:9) GB_parameters(12)];
        elseif length(GB_parameters)==18
            G_parameters=GB_parameters(1:15);       B_parameters=[GB_parameters(1:4) GB_parameters(16) GB_parameters(6:9) GB_parameters(17) GB_parameters(11:14) GB_parameters(18)];
        end%if length(GB_parameters)==6
        G_fit=GB_fit(:,1);                       B_fit=GB_fit(:,2);
        G_residual=GB_residual(:,1);             B_residual=GB_residual(:,2);
        G_l_sq=(G_residual.^2)./(length(freq)-1);%chi-squared calculation
        B_l_sq=(B_residual.^2)./(length(freq)-1);%chi-squared calculation
        [G_parameters,B_parameters]=par_check(G_parameters,B_parameters,freq);
    catch
        G_fit=nan(size(freq,1),size(freq,2));   G_parameters=nan(1,5);  G_l_sq=G_fit;
        B_fit=nan(size(freq,1),size(freq,2));   B_parameters=nan(1,5);  B_l_sq=B_fit;
    end
end%if handles.prefs.simul_peak==0

function [G_parameters,B_parameters]=par_check(G_parameters,B_parameters,freq)
%this function checks to see that the first 5 values in the parmeters
%variable is the most left one, representing the harmonic peak
if length(G_parameters)==10
    if G_parameters(1)-G_parameters(2)<0.9*min(freq)
        G_parameters=[G_parameters(6:10),G_parameters(1:5)];
        B_parameters=[B_parameters(6:10),B_parameters(1:5)];
    elseif G_parameters(1)>G_parameters(6)
        G_parameters=[G_parameters(6:10),G_parameters(1:5)];
        B_parameters=[B_parameters(6:10),B_parameters(1:5)];
    end%if G_parameters(1)>G_parameters(5)
elseif length(G_parameters)==15
    test=[G_parameters(1) G_parameters(6) G_parameters(11)];
    I=find(test==min(test));
    switch I
        case 2
            G_parameters=[G_parameters(6:10) G_parameters(1:5) G_parameters(11:15)];
            B_parameters=[B_parameters(6:10) B_parameters(1:5) B_parameters(11:15)];
        case 3
            G_parameters=[G_parameters(11:15) G_parameters(1:5) G_parameters(6:10)];
            B_parameters=[B_parameters(11:15) B_parameters(1:5) B_parameters(6:10)];
    end%switch I
end%if length(G_parameters)==10

function [fitted_y,residual,parameters]=fit_spectra_con(x0,freq_data,y_data,I,show_GB,lb,ub)%fit spectra to conductance curve
%This function takes the starting guess values ('guess_values'_) and fits a
%Lorentz curve to the the x_data and y_data. The variable 'guess_values' 
%needs to be a 1x5 array. The designation foe each elements is as follows:
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
%Variables, 'lb' and 'ub', defines the lower and upper bound of each of the
%guess_paramters. Both 'lb' and 'ub' are 1x5 array.
if nargin==5
    lb=[0 0 -90 -inf -200];
    ub=[Inf Inf 90 100 200];
end%if nargin==5
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10);
[parameters, ~, ~]=lsqcurvefit(@lfun4c,x0,freq_data(I),y_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
fitted_y=lfun4c(parameters,freq_data);
residual=fitted_y-y_data;
if show_GB==1%checck to see whether or not to show the parameters
    disp('Conductance fitted parameters:');
    disp(parameters');
end%if handles.prefs.show_GB==1

function [fitted_y,residual,parameters]=fit_spectra_sus(x0,freq_data,susceptance_data,I,show_GB,lb,ub)%fit spectra to susceptance curve
%This function takes the starting guess values ('guess_values'_) and fits a
%Lorentz curve to the the x_data and y_data. The variable 'guess_values' 
%needs to be a 1x5 array. The designation foe each elements is as follows:
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
%Variables, 'lb' and 'ub', defines the lower and upper bound of each of the
%guess_paramters. Both 'lb' and 'ub' are 1x5 array.
if nargin==5
    lb=[0 0 -90 -90 -100];
    ub=[Inf Inf 90 100 100];
end%if nargin==5
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10);
[parameters, ~, ~]=lsqcurvefit(@lfun4s,x0,freq_data(I),susceptance_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
fitted_y=lfun4s(parameters,freq_data);
residual=fitted_y-susceptance_data;
if show_GB==1%checck to see whether or not to show the parameters
    disp('Susceptance fitted parameters:');
    disp(parameters');
end%if handles.prefs.show_GB==1

function [fitted_y,residual,parameters]=fit_spectra_both(x0,freq_data,conductance,susceptance,num_peaks,I,handles,lb,ub)
%This function fits both the conductance and susceptance curves simultaneously.
%x0: fitted parameters (see lfun4_both)
%freq_data: frequency array
%conductance: conductance array
%susceptance: susceptance array
%num_peaks: number of peaks to be fitted
%I: indices of conductance and susceptance that will be used for the fitting
%lb: lower bound
%ub: upper bound
if nargin==7
    lb=[.999*min(freq_data) 0 -180 0 -200];
    ub=[1.001*max(freq_data) 2*(max(freq_data)-min(freq_data)) 180 200 200];
end%if nargin==6
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxFunEvals',3e4,'maxiter',3e3);
if length(x0)==6%fitting code for one peak
    [parameters resnorm residual]=lsqcurvefit(@lfun4_both_1,x0,freq_data,[conductance susceptance],[lb,-100],[ub,100],options);
    fitted_y=lfun4_both_1(parameters,freq_data);
    residual=[fitted_y(:,1)-conductance,fitted_y(:,2)-susceptance];
    disp('Fitting 1 peak');
elseif length(x0)==12%fitting code for two peaks
    [parameters resnorm residual]=lsqcurvefit(@lfun4_both_2,x0,freq_data,[smooth(conductance,9) smooth(susceptance,9)],[lb lb -100 -100],[ub ub 100 100],options);
    fitted_y=lfun4_both_2(parameters,freq_data);
    residual=[fitted_y(:,1)-conductance,fitted_y(:,2)-susceptance];
    disp('Fitting 2 peaks');
elseif length(x0)==18%fitting code for three peaks
    [parameters resnorm residual]=lsqcurvefit(@lfun4_both_3,x0,freq_data,[smooth(conductance) smooth(susceptance)],[lb lb lb -100 -100 -100],[ub ub ub 100 100 100],options);
    fitted_y=lfun4_both_3(parameters,freq_data);
    residual=[fitted_y(:,1)-conductance,fitted_y(:,2)-susceptance];
    disp('Fitting 3 peaks');
end% if numpeaks==1

function F_conductance = lfun4c(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_conductance= p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)))+p(5);

function F_susceptance = lfun4s(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_susceptance= -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)))+p(5);

function F_conductance = lfun4c_2(p,x)%this equation fits 2 peaks
F_conductance=lfun4c(p(1:5),x)+lfun4c(p(6:10),x);

function F_susceptance = lfun4s_2(p,x)%this equation fits 2 peaks
F_susceptance=lfun4s(p(1:5),x)+lfun4s(p(6:10),x);

function F_conductance = lfun4c_3(p,x)%this equation fits 3 peaks
F_conductance=lfun4c(p(1:5),x)+lfun4c(p(6:10),x)+lfun4c(p(11:15),x);

function F_susceptance = lfun4s_3(p,x)%this equation fits 3 peaks
F_susceptance=lfun4s(p(1:5),x)+lfun4s(p(6:10),x)+lfun4s(p(11:15),x);

function fcns=lfun4_both_1(p,x)%this function fits 1 peak (cond and sus)
%p(1): f0 maximum frequency (1)                 %p(4): Gmax maximum conductance (1)
%p(2): gamma0 dissipation (1)                     %p(5): Offset value (conductance) (1)
%p(3): phi phase angle difference (1)            %p(6): Offset value (susceptance) (1)
fcns=[lfun4c(p(1:5),x),lfun4s([p(1:4),p(6)],x)];

function fcns=lfun4_both_2(p,x)%this function fits 2 peaks (cond and sus)
%p(1): f0 maximum frequency (1)          p(6): f0 maximum frequency (2)
%p(2): gamma0 dissipation   (1)          p(7): gamma0 dissipation (2)
%p(3): phi phase angle difference (1)    p(8): phi phase angle difference (2)
%p(4): Gmax maximum conductance   (1)    p(9): Gmax maximum conductance (2)
%p(5): Offset value (conductance) (1)    p(10): Offset value (conductance)(2)

%p(11): Offset value (susceptance) (1) p(12): Offset value (susceptance)(2)
fcns=[lfun4c_2(p(1:10),x),lfun4s_2([p(1:4),p(11),p(6:9),p(12)],x)];

function fcns=lfun4_both_3(p,x)%this fucntion fits 3 peaks (cond and sus)
%p(1): f0 maximum frequency (1)          p(6): f0 maximum frequency (2)           p(11): f0 maximum frequency (3)
%p(2): gamma0 dissipation   (1)          p(7): gamma0 dissipation (2)             p(12): gamma0 dissipation (3)
%p(3): phi phase angle difference (1)    p(8): phi phase angle difference(2)      p(13): phi phase angle difference (3)
%p(4): Gmax maximum conductance   (1)    p(9): Gmax maximum conductance (2)       p(14): Gmax maximum conductance (3)
%p(5): Offset value (conductance) (1)    p(10): Offset value (conductance)(2)     p(15): Offset value (conductance (3)

%p(16): Offset value (susceptance) (1)   p(17): Offsetvalue(susceptance)(2)       p(18): offset value (susceptance) (3)
fcns=[lfun4c_3(p(1:15),x),lfun4s_3([p(1:4),p(16),p(6:9),p(17),p(11:14),p(18)],x)];



