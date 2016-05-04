function varargout = two_layer_immersed_v001(varargin)
% Shull research group. Original author: Chyi-Huey Joshua Yeh
% TWO_LAYER_IMMERSED_V001 MATLAB code for two_layer_immersed_v001.fig
%      TWO_LAYER_IMMERSED_V001, by itself, creates a new TWO_LAYER_IMMERSED_V001 or raises the existing
%      singleton*.
%
%      H = TWO_LAYER_IMMERSED_V001 returns the handle to a new TWO_LAYER_IMMERSED_V001 or the handle to
%      the existing singleton*.
%
%      TWO_LAYER_IMMERSED_V001('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWO_LAYER_IMMERSED_V001.M with the given input arguments.
%
%      TWO_LAYER_IMMERSED_V001('Property','Value',...) creates a new TWO_LAYER_IMMERSED_V001 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before two_layer_immersed_v001_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to two_layer_immersed_v001_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help two_layer_immersed_v001

% Last Modified by GUIDE v2.5 25-Dec-2014 20:37:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @two_layer_immersed_v001_OpeningFcn, ...
                   'gui_OutputFcn',  @two_layer_immersed_v001_OutputFcn, ...
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


% --- Executes just before two_layer_immersed_v001 is made visible.
function two_layer_immersed_v001_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to two_layer_immersed_v001 (see VARARGIN)

% Choose default command line output for two_layer_immersed_v001
handles.output = hObject;

%set default values or settings
handles=home_state(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes two_layer_immersed_v001 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles=home_state(handles)
%This function contains all of the default settings of the GUI.
cla(handles.axes1);cla(handles.axes2);
handles.din=[];
%set default values or settings
warning('off','MATLAB:legend:IgnoringExtraEntries');%suppress warning messages associated with havein extra legend entries
handles.din.marker_color=[{[0 0 1]} {[1 0 0]} {[0 0.5 0]} {[0.5 0.5 0]} {[0 0.5 0.5]} {[0.5 0 0.5]}];%set marker colors
handles.din.marker_style=[{['+']} {['o']} {['s']} {['^']} {['v']} {['x']}];%set marker style
handles.din.filepath=pwd;%set default working path directory
handles.din.constants.f1=5e6;%fundamental resonance frequency in Hz
handles.din.constants.zq=8.84e6; %load impedance of quartz, kg/m^2-s
handles.din.constants.del_f_1_e=0.6205;%default error (Hz) 
handles.din.constants.del_f_3_e=1.4309;%default error (Hz) 
handles.din.constants.del_f_5_e=2.4990;%default error (Hz) 
handles.din.constants.del_f_7_e=315 ;%default error (Hz)
handles.din.constants.del_f_9_e=405;%default error (Hz)
handles.din.constants.del_f_11_e=495;%default error (Hz)
handles.din.constants.del_g_1_e=0.3681;%default error (Hz) 
handles.din.constants.del_g_3_e=0.4703;%default error (Hz)
handles.din.constants.del_g_5_e=0.9362;%default error (Hz) 
handles.din.constants.del_g_7_e=50;%default error (Hz)
handles.din.constants.del_g_9_e=50;%default error (Hz)
handles.din.constants.del_g_11_e=50;%default error (Hz)
handles.din.contour.label=0;%default flag for the harm and diss ratios
handles.din.guess_label=1;%default flag for the guess parameters
handles.din.guess_label2=1;%default flag for the guess parameters
handles.din.qcm_label=0;%default plaf for the qcm map
handles.din.contour.res=300;%define resolution for the contour plots
handles.din.contour.harm_ratio_calc=ones(handles.din.contour.res,handles.din.contour.res);%allocate harm_ratio_calc
handles.din.contour.diss_ratio_calc=ones(handles.din.contour.res,handles.din.contour.res);%allocate diss_ratio_calc
handles.din.qcm_map=ones(handles.din.contour.res,handles.din.contour.res);%allocate the matrix for the qcm map
set(handles.harm1,'value',1);
set(handles.harm3,'value',1);
set(handles.radio2,'value',1);
set(handles.d2lam_guess,'string',0.02);%set default value for the d2lam_guess
set(handles.drho_guess,'string',1,'style','text');%set default value for the drho_guess
set(handles.phi_guess,'string',45);%set default value for the phi_guess
set(handles.grho_guess,'string',1e6);%set the default value for grho_guess
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
d2lam_guess_Callback(handles.d2lam_guess, 1, handles);


% --- Outputs from this function are returned to the command line.
function varargout = two_layer_immersed_v001_OutputFcn(hObject, eventdata, handles) 
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
set(handles.status,'string','Status: GUI handles have been reset to its default state!','backgroundcolor','k','foregroundcolor','r');
disp('GUI handles have been reset to its default state!');
handles.din.harmtot=active_harm(handles);%get the active harmonics
try%prompt the user to import the frequency shift values
    set(handles.status,'string','Status: Importing...','backgroundcolor','k','foregroundcolor','r');
    [filename,filepath,~]=uigetfile('.mat','Load frequency shift data',handles.din.filepath);%get the filepath and filename
    handles=import_data(handles,filename,filepath);
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
handles.din.extent=length(timepoints);
set(handles.extent,'string',handles.din.extent);
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


function contour_Callback(hObject, ~, handles)
%This function creates contour plots of the master function.
set(handles.status,'string','Status: Generating contour plots...','backgroundcolor','k','foregroundcolor','r');
disp('Generating contour plots...');
res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
phi_var=linspace(0,90,res);%define range of phi vales
d2lam_var=linspace(0,0.25,res);%define range of d2lam values
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
if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
    try
        drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
        drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
        drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
        drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values
        set(handles.drho_guess,'string',drho_ave*1000);
    catch
        drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
    end
else%if the manual guess for drho has been enabled
    drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
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
contourf(d2lam_var,phi_var,harm_ratio_calc,linspace(-3,3,res),'edgecolor','none');
xlabel(s1,['d/\lambda_',num2str(handles.din.n_ref)],'fontweight','bold','fontsize',14);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',14);
t1=title([num2str(n1),'\Deltaf_',num2str(n2),'/',num2str(n2),'\Deltaf_',num2str(n1)],'fontweight','bold','fontsize',14);
c1=colorbar;
set(c1,'fontsize',12);
set(s1,'fontsize',14);
drawnow;
s2=subplot(1,2,2);%plot the dissipation ratio map
contourf(d2lam_var,phi_var,diss_ratio_calc,linspace(-3,3,res),'edgecolor','none');
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
    sol1=axes;%plot the harmonic ratio solutions    
    plot(sol1,d2lam_dsol,phi_dsol,'.','color','g','linewidth',2);hold(sol1,'on');
    plot(sol1,d2lam_hsol,phi_hsol,'k.','linewidth',2);
    plot(sol1,d2lam_dhsol,phi_dhsol,'m.','linewidth',2);
    sol2=axes;%plot the dissipation ratio solutions
    plot(sol2,d2lam_hsol,phi_hsol,'k.','linewidth',2);hold(sol2,'on');
    plot(sol2,d2lam_dsol,phi_dsol,'.','color','g','linewidth',2);
    plot(sol2,d2lam_dhsol,phi_dhsol,'m.','linewidth',2);    
    %modify the title text for each subplot
    if n1~=1%write out the title text for the harmonic ratio
        set(t1,'interpreter','latex','string',...
            ['\sffamily$\bf\frac{',num2str(n1),'\Delta f_',num2str(n2),'}{',num2str(n2),'\Delta f_',num2str(n1),...
            '}$ = ',num2str(harm_ratio_exp),'$\pm$',num2str(harm_ratio_error)]);
    elseif n1==1
        set(t1,'interpreter','latex','string',...
            ['\sffamily$\bf\frac{\Delta f_',num2str(n2),'}{',num2str(n2),'\Delta f_',num2str(n1),...
            '}$ = ',num2str(harm_ratio_exp),' $\pm$ ',num2str(harm_ratio_error)]);
    end%if n1~=1
    set(t2,'interpreter','latex','string',...
        ['\sffamily$\bf\frac{\Delta \Gamma _',num2str(n3),'}{\Delta f_',num2str(n3),...
        '}$ = ',num2str(diss_ratio_exp),' $\pm$ ',num2str(diss_ratio_error)]);
    pause(1);
    linkaxes([s1,s2, sol1, sol2],'xy');
    set(sol1,'color','none','position',get(s1,'position'),'xlim',get(s1,'xlim'),'ylim',get(s2,'ylim'),'fontsize',14);
    set(sol2,'color','none','position',get(s2,'position'),'xlim',get(s2,'xlim'),'ylim',get(s2,'ylim'),'fontsize',14);
    
end%if get(handles.qcm_map_sol,'value')==1s
drawnow
disp('Plot succesfully generated');



function d2lam_choice_Callback(hObject, eventdata, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
guidata(hObject,handles);

function contour_choice_Callback(hObject, eventdata, handles)
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
guidata(hObject,handles);

function plot_solutions_Callback(hObject, eventdata, handles)
%this function plots out the calculted solutions as a function of time
%NOTE this program currently only plots out data calculated at n1
set(handles.status,'string','Status: Plotting solutions...','backgroundcolor','k','foregroundcolor','r');
disp('Plotting solutions...');

%plot viscoelastic parameters
ff1=figure(4);clf(ff1);
set(ff1,'position',[244,20,1466,414]);
s1=subplot(1,3,1);%plot drho
s2=subplot(1,3,2);%plot grho
s3=subplot(1,3,3);%plot phi
[L_string,xstr]=plot_viscoelastic(s1,s2,s3,handles,1);%plot the viscoelastic properties at the first harmonic (5MHz)
xlabel(s1,xstr,'fontweight','bold','fontsize',14);
ylabel(s1,'d\rho (g/m^2)','fontweight','bold','fontsize',14);
xlabel(s2,xstr,'fontweight','bold','fontsize',14);
ylabel(s2,'|G^*|\rho (Pa-g/cm^3)','fontweight','bold','fontsize',14);
xlabel(s3,xstr,'fontweight','bold','fontsize',14);
ylabel(s3,'\phi (deg.)','fontweight','bold','fontsize',14);
set(s1,'fontsize',14);
set(s2,'fontsize',14);
set(s3,'fontsize',14);
L=legend(s3,'','','','','','');
set(L,'string',L_string);

%plot the raw frequencies
ff2=figure(3);clf(ff2);
set(ff2,'position',[244,564,1466,414]);
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
set(r2,'fontsize',14,'yscale','log');
if get(handles.linkx,'value')==1%if axes linking is turned on
    linkaxes([r1,r2,s1,s2,s3],'x');
end%if get(handles.linkx,'value')==1
set(handles.status,'string','Status: Solutions plotted!','backgroundcolor','k','foregroundcolor','r');
disp('Solutions plotted!');


function qcm_map_Callback(hObject, eventdata, handles)
%this function plots out the contour plots for the master equation
set(handles.status,'string','Status: Generating QCM maps...','backgroundcolor','k','foregroundcolor','r');
disp('Generating QCM maps...');
qcm_map=handles.din.qcm_map;
res=handles.din.contour.res;%number of datapoints in the x and y direction of the plot
phi_var=linspace(0,90,res);%define range of phi vales
d2lam_var=linspace(0,0.5,res);%define range of d2lam values
f1=handles.din.constants.f1;%fundamental resonance freq
zq=handles.din.constants.zq;%quartz  load impedance
f_n=f1*str2double(get(handles.dgliq_harm,'string'));
if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
    try
        drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
        drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
        drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
        drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values
        set(handles.drho_guess,'string',drho_ave*1000);
    catch
        drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
    end
else%if the manual guess for drho has been enabled
    drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
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
contourf(d2lam_var,phi_var,real(qcm_map),linspace(-3,3,300),'edgecolor','none');
xlabel(s1,['d/\lambda_',num2str(n)],'fontweight','bold','fontsize',12);
ylabel(s1,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Deltaf_',num2str(n),'/\Deltaf_s_',num2str(n)],'fontweight','bold','fontsize',12);
colormap(jet); colorbar;
waitbar(.95,h);
s2=subplot(1,2,2);%plot the imaginary component
contourf(d2lam_var,phi_var,imag(qcm_map),linspace(-3,3,300),'edgecolor','none');
xlabel(s2,['d/\lambda_',num2str(n)],'fontweight','bold','fontsize',12);
ylabel(s2,'\phi (deg.)','fontweight','bold','fontsize',12);
title(['\Delta\Gamma_',num2str(n),'/\Deltaf_s_',num2str(n)],'fontweight','bold','fontsize',12);
colormap(jet); colorbar;
waitbar(1,h);
cf=z_star_liq_n/(drho_ave*f_n);
set(ff,'name',['drho=',num2str(drho_ave),'g/m^2  cf=',num2str(cf)]);
set(handles.status,'string','Status: Plot succesfully generated!','foregroundcolor','r','backgroundcolor','k');
try  delete(h); end%try
disp('Plot succesfully generated');

function [L_string,xstr]=plot_viscoelastic(s1,s2,s3,handles,n1)
%this function plots the viscoelastic parrameters in subplots s1, s2, and
%s3 at harmonic n1
sols=fieldnames(handles.din.stored_solutions);%extract out the filenames (n_<harmonic at which the parameters were calcaulated at>_<harmonics used to calc data>
count=1;
L_string=[];
for dum=1:length(sols)
    index=str2double(sols{dum}(3));
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
    if get(handles.norm_time,'value')==1%run this cord if the option to normalize the timepoints is turned on
        xdata=xdata-min(xdata);
    end%if get(handles.norm_time,'value')==1
    L_string=[L_string,{temp(end-2:end)}];%string for legend box
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
    if get(handles.norm_time,'value')==1%run this cord if the option to normalize the timepoints is turned on
        xdata=xdata-xdata(1);
    end%  if get(handles.norm_time,'value')==1
    plot(r1,xdata(:,1),freq_shifts(:,handles.din.harmtot(dum)+1),'-','color',handles.din.marker_color{dum},'linewidth',1);hold(r1,'on');%plot experimental Df
    plot(r2,xdata(:,1),freq_shifts(:,handles.din.harmtot(dum)+2),'-','color',handles.din.marker_color{dum},'linewidth',1);hold(r2,'on');%plot pred. Df
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
        index=str2double(sols{dum}(3));
            if index==dum1
                index2(count)=dum;
                count=count+1;
            end%if index==1
    end%for dum=1:length(sols)
    for dum=1:length(index2)%calc for each dataset combination
        temp=sols{index2(dum)};
        [xstr,xdata]=determine_time_units(handles.din.stored_solutions.(temp).g_shifts(:,1),handles);%determine xdata
        if get(handles.norm_time,'value')==1%run this cord if the option to normalize the timepoints is turned on
            xdata=xdata-min(xdata);
        end%  if get(handles.norm_time,'value')==1
        ydata1=handles.din.stored_solutions.(temp).f_shifts(:,3);%calc Df
        ydata2=handles.din.stored_solutions.(temp).g_shifts(:,3);%Calc Dg       
        if dum1==1
            L_string=[L_string,{[temp(end-2:end)]}];%string for legend box
        end%if dum1==1
        if flag==0%if the flag is turned off, plot the pred freq black
            hold(r1,'on');hold(r1,'on');
            plot(r1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none');
            plot(r2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color','k','linestyle','none'); 
        else
            index3=0.5*(str2double(sols{index2(dum)}(3))+1);
            hold(r1,'on');hold(r1,'on');
            plot(r1,xdata,ydata1,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
            plot(r2,xdata,ydata2,'linewidth',2,'markersize',14,'marker',marker_style{dum},'color',marker_color{index3},'linestyle','none');
        end%if flag==0
    end%for dum=1:length(index2)
end%for dum1=1:2:11



function qcm_map_sol_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION CONTAINS CALLBACK FUNCTIONS IN THE ACTIVE HARMONICS PANEL

function harm1_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
function harm3_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
function harm5_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
function harm7_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
function harm9_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);
function harm11_Callback(hObject, eventdata, handles)
handles.din.harmtot=active_harm(handles);%get the active harmonics
handles=import_data(handles,handles.din.filename,handles.din.filepath);%reload the data
handles=plot_fcn(handles,handles.din.rawdata);%plot the rawdata
guidata(hObject,handles);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS SECTION DEALS WITH CALLBACK FUNCTIONS ASSOCIATED WITH THE CALCULATION
%OPTIONS PANEL

function radio1_Callback(hObject, eventdata, handles)
function radio2_Callback(hObject, eventdata, handles)
function radio3_Callback(hObject, eventdata, handles)
function radio4_Callback(hObject, eventdata, handles)
function radio5_Callback(hObject, eventdata, handles)
function radio6_Callback(hObject, eventdata, handles)
function harmchoice1_Callback(hObject, eventdata, handles)
function harmchoice2_Callback(hObject, eventdata, handles)
function harmchoice3_Callback(hObject, eventdata, handles)
function harmchoice4_Callback(hObject, eventdata, handles)
function harmchoice5_Callback(hObject, eventdata, handles)
function harmchoice6_Callback(hObject, eventdata, handles)

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
guidata(hObject,handles);

function edit_drho_Callback(hObject, ~, handles)
%This function allows for the user to manually guess the value for drho
if get(handles.edit_drho,'value')==1%if the value equals 1, allow for the guess value of drho to be editable
    set(handles.drho_guess,'style','edit');
else
    set(handles.drho_guess,'style','text');
end%if get(handles.edit_drho,'value')==1


function cursor_Callback(hObject, ~, handles)
%This function allows the user to pick a point from the plots. The
%datapoint is then used for analysis.
handles.din.stored_solutions=[];%clear out the stored solutions
button=1;%this is a flag that represent the which button was pushed during datapoint selection
count=1;
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
            try%sometimes interp1 does not always work
                handles.din.cursor.(namefi)=interp1(xdata_f,ydata_f,selected_time0,'linear');%interpolate the frequency shift defined by selected_time0
                handles.din.cursor.(namegi)=interp1(xdata_g,ydata_g,selected_time0,'linear');%interpolate the frequency shift defined by selected_time0
            catch%if it does not work try a different method (find nearest timepoint);
                ind3=find(xdata_f-selected_time0==min(abs(xdata_f-selected_time0)),1,'first');%find the index associated with the nearest timepoint
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

function solve_all_Callback(hObject, eventdata, handles)
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
        waitbar(dum/extent,h);
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
                handles=qcm_solver(handles,count);%solve for viscoelastic parameters  
                count=count+1;
            end%if isnan(handles.din.cursor.(namefi))==1||isnan(handles.din.cursor.(namegi))==1      
        catch err_msg
            assignin('base','err_msg',err_msg);
            set(handles.status,'string','Status: Unable to find solution!','backgroundcolor','k','foregroundcolor','r');
            disp('Unable to find solution!');
        end%try
    end%dum2=1:length(handles.din.harmtot)    
end%for dum=length(timepoints)
try delete(h); end%try to delete the status bar
guidata(hObject,handles);
set(handles.status,'string','Status: Calculation finished!','backgroundcolor','k','foregroundcolor','r');
disp('Calculation finished!');

function solve_inc_Callback(hObject, eventdata, handles)


function auto_solve_Callback(hObject, eventdata, handles)


function drho_guess_Callback(hObject, eventdata, handles)
%This function recalculates the 3rd guess parameter (that is in red).
grho_state=get(handles.grho_guess,'foregroundcolor');%get the color of the font
d2lam_state=get(handles.d2lam_guess,'foregroundcolor');%get the color of the font
if sum(grho_state)==1&&sum(d2lam_state)==0%recalculate grho_guess if grho_guess is red
    handles=d2lam_guess_Callback(handles.d2lam_guess, eventdata, handles);
elseif sum(grho_state)==0&&sum(d2lam_state)==1%recalculate d2lam_guess if d2lam_giess is red
    handles=grho_guess_Callback(handles.grho_guess, eventdata, handles);
end%if sum(grho_state)==1&&sum(d2lam_state)==0
handles.din.guess_label=rand(1);%"relabel" the guess parameters, this signals the contour callback function to recalculate the harm/diss ratio matrices
handles.din.guess_label2=rand(1);
guidata(hObject,handles);


function phi_guess_Callback(hObject, eventdata, handles)
%This function recalculates the 3rd guess parameter (that is in red).
grho_state=get(handles.grho_guess,'foregroundcolor');%get the color of the font
d2lam_state=get(handles.d2lam_guess,'foregroundcolor');%get the color of the font
if sum(grho_state)==1&&sum(d2lam_state)==0%recalculate grho_guess if grho_guess is red
    handles=d2lam_guess_Callback(handles.d2lam_guess, eventdata, handles);
elseif sum(grho_state)==0&&sum(d2lam_state)==1%recalculate d2lam_guess if d2lam_guess is red
    handles=grho_guess_Callback(handles.grho_guess, eventdata, handles);
end%if sum(grho_state)==1&&sum(d2lam_state)==0
guidata(hObject,handles);


function handles=grho_guess_Callback(hObject, ~, handles)
%Running this callback function will calculate the predicted d2lambda value.
phi=str2double(get(handles.phi_guess,'string'));%extract out the value for the phi_guess
drho=str2double(get(handles.drho_guess,'string'))./1000;%extract out the value for the drho_guess
grho=str2double(get(handles.grho_guess,'string')).*1000;%extract out the value for the grho_guess
f1=handles.din.constants.f1;%fundamental resonance frequency
d2lam=(drho.*f1.*1.*cosd(phi./2))./(sqrt(grho));%calc. the d2lam value
set(handles.d2lam_guess,'string',d2lam,'foregroundcolor','r');
set(handles.grho_guess,'foregroundcolor','k');
guidata(hObject,handles);


function handles=d2lam_guess_Callback(hObject, ~, handles)
%Running this callback function will calculate the predicted Grho value.
d2lam=str2double(get(handles.d2lam_guess,'string'));%extract out the value for the d2lam_guess
phi=str2double(get(handles.phi_guess,'string'));%extract out the value for the phi_guess
drho=str2double(get(handles.drho_guess,'string'))./1000;%extract out the value for the drho_guess
f1=handles.din.constants.f1;%fundamental resonance frequency
grho=(((drho./d2lam).*f1.*1.*cosd(phi./2)).^2)./1000;%calc. the grho value
set(handles.grho_guess,'string',grho,'foregroundcolor','r');
set(handles.d2lam_guess,'foregroundcolor','k');
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

function harm_ratio_error=h_ratio_error(del_f_n2,del_f_n1,del_f_n1_e, del_f_n2_e,n1,n2)
%this function calculates the error assosciated with the harmonic ratio.
%This error calculation is based on standard propagation of error. See
%DeNolf et al. Langmuir, 2014 for more details on error analysis.
harm_ratio_error=sqrt((((n1.*del_f_n2.*del_f_n1_e)./(n2.*(del_f_n1.^2))).^2)+(((n1.*del_f_n2_e)./(n2.*del_f_n1)).^2));

function diss_ratio_error=d_ratio_error(del_g_n3,del_f_n3,del_g_n3_e,del_f_n3_e)
%This function calcualtes the error associated with the dissipation ratio.
%This error calculation is based on standard propgation of error. See
%DeNolf et al. Langmuir, 2014 for more details on error analysis.
diss_ratio_error=sqrt((((del_g_n3_e)./((del_f_n3))).^2)+(((del_g_n3.*del_f_n3_e)./(del_f_n3.^2)).^2));

function [harm_ratio_error,diss_ratio_error,handles]=calc_ratio_error(handles,n1,n2,n3,del_f_n2,del_f_n1,del_g_n3,del_f_n3)
%This function calculates the error for the harmonic and dissipation error
%based in the experimental error from the frequency shifts.
%handles: gui handles structure
%n1: first harmonic dataset used to calculate the harmonic ratio
%n2: second harmonic dataset used to calculate the harmonic ratio
%n3:  third harmonic dataset used to calcualte the dissipation ratio
%calculate error in ratios
del_f_n1_e=handles.din.constants.(['del_f_',num2str(n1),'_e']);
del_f_n2_e=handles.din.constants.(['del_f_',num2str(n2),'_e']);
del_g_n3_e=handles.din.constants.(['del_g_',num2str(n3),'_e']);
del_f_n3_e=handles.din.constants.(['del_f_',num2str(n3),'_e']);
harm_ratio_error=h_ratio_error(del_f_n2,del_f_n1,del_f_n1_e, del_f_n2_e,n1,n2);%error in the harmonic ratio
diss_ratio_error=d_ratio_error(del_g_n3,del_f_n3,del_g_n3_e,del_f_n3_e);%error in the dissipation ratio


function F=immersed_2layer_solver(handles,n1,n2,n3,n,z_star_liq_n,harm_ratio_exp,diss_ratio_exp,guess_values)
%This function is the function that is called by fsolve. It calculates the
%dissipation and harmonic ratio and takes the difference relative to the
%experimental values (F).
%n1: harmonic at n1
%n2: harmonic at n2
%n: hamonic at n or the harmonic in which the liq load impedance is
%associated with
%z_star_liq_n: liq load impedance at harmonic n
%harm_ratio_exp: experimentally obtained harmonic ratio
%diss_ratio_exp: experimentally obtained dissipation ratio
%guess_values(1): d2lam_n1 guess
%guess_values(2): phi
f1=handles.din.constants.f1;
zq=handles.din.constants.zq;
d2lam_n1=guess_values(1);%guess d2lambda ratio at harmonic n1
phi=guess_values(2);%guess viscoelastic phase angle
d2lam_n2=d2lam_n2_calc(d2lam_n1,n1,n2,phi);%based on guess, calc. d2lambda at harmonic n2
d2lam_n3=d2lam_n2_calc(d2lam_n1,n1,n3,phi);%based on guess, calc. d2lambda at harmonic n3
z_star_liq_n1=zqliq_n2_calc(z_star_liq_n,n,n1);%liq load impedance at harmonic n1
z_star_liq_n2=zqliq_n2_calc(z_star_liq_n,n,n2);%liq load impedance at harmonic n2
z_star_liq_n3=zqliq_n2_calc(z_star_liq_n,n,n3);%liq load impedance at harmonic n3
if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
    try
        drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
        drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
        drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
        drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values
        if isnan(drho_ave)==0
            set(handles.drho_guess,'string',drho_ave*1000);
        end%if isnan(drho_ave)==0
    catch
        drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
    end
else%if the manual guess for drho has been enabled
    drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
end%if get(handles.edit_drho,'value')==0
f_n1=f1*n1;%resonant frequency at n1 in Hz
f_n2=f1*n2;%resonant frequency at n2 in Hz
f_n3=f1*n3;%resonant frequency at n3 in Hz
%only used drho_n2 because drho_n1, n2, n3 should be similar
harm_ratio_calc=real(master(d2lam_n2,phi,z_star_liq_n2,f_n2,drho_ave))./real(master(d2lam_n1,phi,z_star_liq_n1,f_n1,drho_ave));%calculate the harmonic ratio
diss_ratio_calc=imag(master(d2lam_n3,phi,z_star_liq_n3,f_n3,drho_ave))./real(master(d2lam_n3,phi,z_star_liq_n3,f_n3,drho_ave));%calculate the dissipation ratio
% drho_calc=100.*(drho_recalc(f_n1,zq,f1,n1,master(d2lam_n1,phi,z_star_liq_n1,f_n1,drho_ave))-drho_ave);
F=[harm_ratio_calc-harm_ratio_exp;diss_ratio_calc-diss_ratio_exp];%calculates the difference between the calculated and experimental ratios

function handles=qcm_solver(handles,count)
%this function calculates the two variables, d2lambda and phi.
%handles: the gui handles structure
%count: keeps track which row to store the calculated viscoelasetic
%parameters in a variable
[radiotot,harmchoice]=radio_total(handles);%determine how many datasets and which harmonics that will be used to perform the calculation
f1=handles.din.constants.f1;%fundamentalr resonance freq in Hz
zq=handles.din.constants.zq;% the quartz load impedance
for dum=1:length(radiotot)
    calc_table=nan(6,8);    %preallocate the calculation table
    handles.din.solved.harmchoice=harmchoice(dum);
    [n1,n2,n3]=determine_harm_choice(harmchoice(dum));%determine what harmonic datasets will be used to calculate viscoelastic parameters
    %calculate the liquid load impedance
    dgliq_n=str2double(get(handles.dgliq,'string'));%diss shift used for liq load impedance calc
    dfliq_n=-dgliq_n;%freq shift used for liq load impedance calc
    z_star_liq_n=zqliq_n_calc(dfliq_n,dgliq_n,f1,zq);%calc the liq load impedance at dgliq_harm
    dgliq_harm=str2double(get(handles.dgliq_harm,'string'));%harmonic in which the liq load impedance was calc at
    %extract out the guess parameters
    guess_values(1)=str2double(get(handles.d2lam_guess,'string'));%guess value for d2lambda at n1
    guess_values(2)=str2double(get(handles.phi_guess,'string'));%guess value for phi
   %calculate experimental harmonic and dissipation ratio
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
    try
        options=optimset('Display','off','tolfun',1e-6,'tolx',1e-6,'Maxiter',100000,'MaxFunEvals',100000,'Algorithm','levenberg-marquardt');
        [solved,fval,exitflag,output,jacobian]=fsolve(fcns,guess_values,options);
        output.fval=fval;%store the final difference betwen the harm and diss ratios from exp and harm and diss rations from fsolve
        output.jacobian=jacobian;%store the jacobian matrix from fsolve
        output.exitflag=exitflag;%store the exitflag that was outputed from fsolve
        disp(['exitflag: ',num2str(exitflag)]);
        if exitflag~=1%stop the code if the exitflag is not 1
            set(handles.status,'string','Status: Exitflag is not 1!','backgroundcolor','y','foregroundcolor','r');
            disp('Exitflag is not 1!');
            return
        else
            set(handles.status,'string','Calculating...','backgroundcolor','k','foregroundcolor','r');
        end%if exitflag~=1
        handles.din.solved.output_misc_from_fsolve=output;%store the output structure in the handles structure
        handles.din.solved.solution=solved;%store the solution, solved(1): d2lam at n1 and solved(2): phi
        %solve for the other viscoelastic parameters        
        if get(handles.edit_drho,'value')==0% if the manual guess for drho has been disabled
            try
                drho_n1=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n1)]),n1,f1,zq);%areal mass calc. based on harmonic n1 (not used)
                drho_n2=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n2)]),n2,f1,zq);%areal mass calc. based on harmonic n2
                drho_n3=drho_est(handles.din.cursor.(['interp_harmfi',num2str(n3)]),n3,f1,zq);%areal mass calc. based on harmonic n3 (not used)
                drho_ave=mean(unique([drho_n1,drho_n2,drho_n3]));%take the average value from all of the drho values
                set(handles.drho_guess,'string',drho_ave);
            catch
                drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%kg/m^2
            end
        else%if the manual guess for drho has been enabled
            drho_ave=str2double(get(handles.drho_guess,'string'))./1000;%in kg/m^2
        end%if get(handles.edit_drho,'value')==0
            unique_n=[n1;n2;n3;active_harm(handles)];
            unique_n=unique(unique_n,'stable');
        for dum=1:length(unique_n)
            name=['n_',num2str(unique_n(dum)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
            handles=calc_viscoelastic_par(handles,unique_n(dum),f1,z_star_liq_n,zq,solved,drho_ave,dgliq_harm,name,n1);%calculate the viscoelastic parameters
            
            %store the calculated values in a matrix that will be used to output the results in the handles.uitable2
            calc_table((unique_n(dum)+1)/2,:)=[handles.din.solved.(name).drho,handles.din.solved.(name).grho,handles.din.solved.(name).phi,...
                handles.din.solved.(name).d2lam,handles.din.solved.(name).lam_rho,real(handles.din.solved.(name).norm_delfstar),...
                imag(handles.din.solved.(name).norm_delfstar),abs(handles.din.solved.(name).cf)];
            %store the predicted freq shifts in a matrix that will be used to output the results in the handles.uitable3
            pred_shifts((unique_n(dum)+1)/2,:)=[handles.din.cursor.(['interp_harmfi',num2str(unique_n(dum))]),real(handles.din.solved.(name).complex_pred_freq),...
                handles.din.cursor.(['interp_harmgi',num2str(unique_n(dum))]),imag(handles.din.solved.(name).complex_pred_freq)];
            
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
    end%try
    try%calculate the error
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp+harm_ratio_error,diss_ratio_exp+diss_ratio_error,x);
        [solved_e1,~,exitflage1,~,~]=fsolve(fcns,guess_values,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp-harm_ratio_error,diss_ratio_exp-diss_ratio_error,x);
        [solved_e2,~,exitflage2,~,~]=fsolve(fcns,guess_values,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp+harm_ratio_error,diss_ratio_exp-diss_ratio_error,x);
        [solved_e3,~,exitflage3,~,~]=fsolve(fcns,guess_values,options);
        fcns=@(x)immersed_2layer_solver(handles,n1,n2,n3,dgliq_harm,z_star_liq_n,harm_ratio_exp-harm_ratio_error,diss_ratio_exp+diss_ratio_error,x);
        [solved_e4,~,exitflage4,~,~]=fsolve(fcns,guess_values,options);
        error0.solved_e1=solved_e1;
        error0.solved_e2=solved_e2;
        error0.solved_e3=solved_e3;
        error0.solved_e4=solved_e4;
        for dum0=1:4
            for dum3=1:length(unique_n)
                name=['n_',num2str(unique_n(dum3)),'_',num2str(n1),num2str(n2),num2str(n3)];%nomenclature: n_<harmonic in which values are calc. at>_<harm used to calc properties>
                name_e=['solved_e',num2str(dum0)];
                handles=calc_viscoelastic_par(handles,unique_n(dum3),f1,z_star_liq_n,zq,error0.(name_e),drho_ave,dgliq_harm,name_e,n1);%calc the viscoelastic parameters
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
        set(handles.d2lam_guess,'string',handles.din.solved.(name).d2lam);
        set(handles.drho_guess,'string',handles.din.solved.(name).drho);
        set(handles.phi_guess,'string',handles.din.solved.(name).phi);
        set(handles.grho_guess,'string',handles.din.solved.(name).grho);%in units of Pa-g/cm^2
    end%get(handles.static_guess)==0
end%for dum=1:length(radiotot)

function handles=calc_viscoelastic_par(handles,n1,f1,z_star_liq_n,zq,solved,drho_ave,dgliq_harm,name,n0)
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
handles.din.solved.(name).norm_delfstar=master(handles.din.solved.(name).d2lam,solved(2),zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1),f1*n1,drho_ave);%normalized complex frequency shift at harmonic n1
if sum([nc1==n1,nc2==n1,nc3==n1])==0%if n1 is not one of the harmonics that are being used to calculate the viscoeasltic parameters, use the calculated areal mass at nc1 to predict the frequency shift. 
    %This prevents the predicted freq shift from being the same as the exp. freq shift    
    name0=['n_',num2str(nc1),'_',num2str(nc1),num2str(nc2),num2str(nc3)];%nomenclature: n_<n1>_<harm used to calc properties>
    handles.din.solved.(name).drho=handles.din.solved.(name0).drho;%calculate drho for harmonic nc1 that was used in g/m^2
    handles.din.solved.(name).complex_pred_freq=handles.din.solved.(name).norm_delfstar.*sauerbrey(n1,f1,handles.din.solved.(name).drho./1000,zq);%calculate the complex frequency shifts
else %
    handles.din.solved.(name).drho=drho_calc(handles.din.cursor.(['interp_harmfi',num2str(n1)]),zq,f1,n1,handles.din.solved.(name).norm_delfstar).*1000;%calculate drho for each harmonic that was used in g/m^2
    handles.din.solved.(name).complex_pred_freq=handles.din.solved.(name).norm_delfstar.*sauerbrey(n1,f1,handles.din.solved.(name).drho./1000,zq);%calculate the complex frequency shifts
end%if sum([nc1==n1,nc2==n1,nc3==n1])==0
handles.din.solved.(name).phi=solved(2);%store the viscoelastic phase angle at n1 (degrees)
handles.din.solved.(name).grho=grho_calc(handles.din.solved.(name).drho./1000,handles.din.solved.(name).d2lam,n1*f1,solved(2))./1000;%calculate grho at n1 in units of Pa-g/cm^3
handles.din.solved.(name).cf=zqliq_n2_calc(z_star_liq_n,dgliq_harm,n1)./(n1.*f1.*handles.din.solved.(name).drho);%calculate the "cf" term
handles.din.solved.(name).lam_rho=sqrt(handles.din.solved.(name).grho)./(n1.*f1.*cosd(solved(2)/2));%calculate \lambda\rho


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
    end%switch get(handles.time_units,'value')


% --------------------------------------------------------------------
function paintbrush_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to paintbrush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
brush
