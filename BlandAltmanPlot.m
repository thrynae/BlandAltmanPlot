function h=BlandAltmanPlot(var1, var2, varargin)
%Create a Bland-Altman plot and return a struct with results and handles.
%This function does not require any toolbox.
%
%Syntax:
% BlandAltmanPlot()
%    Generate an example in the current axes.
% BlandAltmanPlot(var1,var2)
%    Create a Bland-Altman plot with the two entered variables. All other settings are set to the
%    default, as described below. Both must be numeric vectors with 2 dimensions or fewer (so
%    1x1x30 sized arrays will return an error). The y-values are var2-var1, the x-values are
%    (var1+var2)/2, unless the plot_x_mean switch is set to false, in which case var1 determines
%    the x-values.
% BlandAltmanPlot(___,Name,Value)
%    Change optional parameters. See below for the list of options that can be changed.
% BlandAltmanPlot(___,optionstruct)
%    Change optional parameters. Each field should be one of the parameters detailed below. Missing
%    fields are filled with the default.
% h=BlandAltmanPlot(___)
%    Returns a struct with values and object handles. This can be used to extract calculated
%    values, to delete specific objects, or to change the plot type.
%
%Parameters:
%AddDetailsText      [default=true] Add the details to the plot as text. This will add the mean and
%                    limits of agreement above the plotted lines. The confidence intervals will be
%                    added under the plotted lines if enabled by the plotCI parameter. The text
%                    elements are placed on the far right side of the plot.
%alpha               [default=0.05] The alpha value is used for the limits of agreement, as well as
%                    for the confidence intervals.
%plotCI              [default=true] Add the CIs of the mean and LoAs to the plot with the errorbar
%                    function. If AddDetailsText is set to true, this parameter controls both the
%                    whiskers and the text. The whiskers are plotted close to the left side of the
%                    plot.
%GroupIndex          [default=[]] You can split the plot into an array of line objects so you can
%                    adjust properties like Color and Marker for subgroups. Setting this to an
%                    empty array will keep all points in a single group.
%                    This parameter must contain only positive integer values, as it is used to
%                    index h.plot.data.
%plot_x_mean         [default=true] If set to false, the first input is used for the x-coordinate,
%                    instead of using the mean.
%StoreToAppdata      [default=true] If set to true, the output of this function will be stored in
%                    the axes object with setappdata, so it can be easily retrieved. The output
%                    struct is set to the data field named 'HJW___BlandAltmanPlot___data'.
%Parent              [default=[]] This determines the parent axes for the plots and texts. If left
%                    empty, gca is used to determine the target axes. If TargetAxes and Parent both
%                    are non-empty, TargetAxes is ignored.
%TargetAxes          [default=[]] This is the same as the Parent switch. If both are non-empty,
%                    TargetAxes is ignored. This switch is maintained for backwards compatibility.
%TextDigitsDisplayed [default=4] This is either the number of digits used in the text elements, or
%                    a FormatSpec that makes num2str return a non-empty char.
%xxyy                [default=[]] This parameter controls the axis range. It has to be a 4-element
%                    vector. Any NaN values will be auto-determined based on the input data (taking
%                    into account the value of plot_x_mean). An empty input is equivalent to [NaN
%                    NaN NaN NaN].
%
%Output:
%h.input.alpha       - alpha value used for plots and calculations
%       .var1        - first variable used for plot
%       .var2        - second variable used for plot
% .data.mu           - mean value of var2-var1
%      .loa          - lower and upper limits of agreement
%      .CI.mu        - lower and upper bound of the CI of the mean
%         .loa_lower - lower and upper bound of the CI of the lower LoA
%         .loa_upper - lower and upper bound of the CI of the upper LoA
% .xxyy              - axis extents used
% .plot.data         - handle to the data plot object (may be an array)
%      .mean         - handle to the mean line
%      .loa_lo       - handle to the lower LoA line
%      .loa_hi       - handle to the upper LoA line
%      .CI_lo        - handle to the lower LoA errorbar
%      .CI_mu        - handle to the mean errorbar
%      .CI_hi        - handle to the upper LoA errorbar
% .text.mean         - handle to the 'mean: %f' text object
%      .meanCI       - handle to the 'CI: %f-%f' text object
%      .loa_hi       - handle to the 'LoA: %f' text object
%      .loa_hi_CI    - handle to the 'CI: %f-%f' text object
%      .loa_lo       - handle to the 'LoA: %f' text object
%      .loa_lo_CI    - handle to the 'CI: %f-%f' text object
%
%The method for calculating characteristic values and the example data were taken from
%Bland&Altman (Lancet, 1986, i:307-310) http://dx.doi.org/10.1016/S0140-6736(86)90837-8
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.2.1                                                         |%
%|  Date:    2021-06-11                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - The implementation of tinv and alpha_to_Z are by Star Strider. These implementations mean you
%   can use this function without the Statistics Toolbox.
%   https://www.mathworks.com/matlabcentral/fileexchange/56500
%   https://www.mathworks.com/matlabcentral/answers/45173#answer_55318

if nargin<2
    if nargin==0
        %example data taken from Bland&Altman (Lancet, 1986, i:307-310)
        % http://dx.doi.org/10.1016/S0140-6736(86)90837-8
        %(only first measurement of each peak flow meter)
        var1=[494,395,516,434,476,557,413,442,650,433,417,656,267,478,178,423,427];
        var2=[512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451];
    else
        error('Incorrect number of input argument.')
    end
end
[success,options,ME]=BlandAltmanPlot_parse_inputs(var1,var2,varargin{:});
if ~success
    %The throwAsCaller function was introduced in R2007b, hence the rethrow here.
    rethrow(ME)
else
    [StoreToAppdata, alpha, plot_x_mean, plotCI, AddDetailsText, ...
        h_ax, TextDigitsDisplayed, xxyy, GroupIndex]=...
        deal(options.StoreToAppdata, options.alpha, options.plot_x_mean, options.plotCI, ...
        options.AddDetailsText, options.TargetAxes, options.TextDigitsDisplayed, options.xxyy, ...
        options.GroupIndex);
end

%Compute Bland-Altman plot characteristics.
[mu,loa,CI]=BlandAltman_values(var1,var2,alpha);

h=struct;%Store all handles to this struct.

%Prepare first part of output.
h.input.alpha=alpha;
h.input.var1=var1;
h.input.var2=var2;
h.data.mu=mu;
h.data.loa=loa;
h.data.CI=CI;

%Determine x and y data based on plot_x_mean.
if ~isequal(size(var1),size(var2))
    %Reshape to vectors.
    var1=var1(:);var2=var2(:);
end
if plot_x_mean
    x=(var1+var2)/2;
    y=var2-var1;
else
    x=var1;
    y=var2-var1;
end

%Store the value of xxyy in the output struct, so it can be used to round the axis extent to
%appropriate values on a case by case basis.
h.xxyy=[min(x) max(x) min(y) max(y)];
%Consider limits of agreement CI when calculating the bounds.
h.xxyy(3)=min(h.xxyy(3),CI.loa_lower(1));
h.xxyy(4)=max(h.xxyy(4),CI.loa_upper(2));
%Add a small margin so the points are not on the exact edge of the plot.
margin=5;%in percent
meanxxyy=[[1 1]*mean(h.xxyy(1:2)) [1 1]*mean(h.xxyy(3:4))];
h.xxyy=meanxxyy+(h.xxyy-meanxxyy)*(1+margin/100);
if isempty(xxyy)
    xxyy=[NaN NaN NaN NaN];
end
%Replace any NaN by the auto-determined extent.
xxyy(isnan(xxyy))=h.xxyy(isnan(xxyy));h.xxyy=xxyy;

%Plot the data, mean and LoAs.
NextPlot=get(h_ax,'NextPlot');%Get property to retain state.
set(h_ax,'NextPlot','add')
if ~isempty(GroupIndex)
    Color=[];
    for ind=unique(GroupIndex(:)).'
        L=GroupIndex==ind;
        h.plot.data(ind)=plot(x(L),y(L),'.','Parent',h_ax);
        if isempty(Color)
            Color=get(h.plot.data(ind),'Color');
        else
            set(h.plot.data(ind),'Color',Color);
        end
    end
else
    h.plot.data=plot(x,y,'.','Parent',h_ax);
end
h.plot.mean=plot(xxyy(1:2),mu*[1 1],'k','Parent',h_ax);
h.plot.loa_lo=plot(xxyy(1:2),loa(1)*[1 1],'k--','Parent',h_ax);
h.plot.loa_hi=plot(xxyy(1:2),loa(2)*[1 1],'k--','Parent',h_ax);
if plotCI
    err1=diff(CI.loa_lower)/2;
    err2=diff(CI.mu)/2;
    xe=xxyy(1)+0.05*diff(xxyy(1:2));%Shift slightly right.
    h.plot.CI_lo=errorbar(xe,loa(1),err1,'ko');
    h.plot.CI_mu=errorbar(xe,mu    ,err2,'ko');
    h.plot.CI_hi=errorbar(xe,loa(2),err1,'ko');
end
set(h_ax,'NextPlot',NextPlot)%Restore state (probably hold('off') ).
axis(h_ax,xxyy)
if AddDetailsText
    %Add text elements.
    %Use num2str to make use of automatic choice of FormatSpec for a fixed number of digits. This
    %is guaranteed to return a non-empty char array, as that was tested in the parser.
    d=TextDigitsDisplayed;%Use a shorter variable name for better readability of the code.
    
    %Add mean text.
    tx=xxyy(2);ty=mu;
    txt=sprintf('mean: %s',num2str(ty,d));
    h.text.mean=text(tx,ty,txt,'Parent',h_ax,...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    if plotCI
        txt=sprintf('CI: %s - %s',...
            num2str(CI.mu(1),d),num2str(CI.mu(2),d));
        h.text.meanCI=text(tx,ty,txt,'Parent',h_ax,...
            'HorizontalAlignment','right','VerticalAlignment','top');
    end
    
    %Add upper limit of agreement.
    ty=loa(2);
    txt=sprintf('LoA: %s',num2str(ty,d));
    h.text.loa_hi=text(tx,ty,txt,'Parent',h_ax,...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    if plotCI
        txt=sprintf('CI: %s - %s',...
            num2str(CI.loa_upper(1),d),num2str(CI.loa_upper(2),d));
        h.text.loa_hi_CI=text(tx,ty,txt,'Parent',h_ax,...
            'HorizontalAlignment','right','VerticalAlignment','top');
    end
    
    %Add lower limit of agreement.
    ty=loa(1);
    txt=sprintf('LoA: %s',num2str(ty,d));
    h.text.loa_lo=text(tx,ty,txt,'Parent',h_ax,...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
    if plotCI
        txt=sprintf('CI: %s - %s',...
            num2str(CI.loa_lower(1),d),num2str(CI.loa_lower(2),d));
        h.text.loa_lo_CI=text(tx,ty,txt,'Parent',h_ax,...
            'HorizontalAlignment','right','VerticalAlignment','top');
    end
end

if StoreToAppdata
    setappdata(h_ax,'HJW___BlandAltmanPlot___data',h)
end
if nargout==0
    clear h
end
end
function [mu,loa,CI]=BlandAltman_values(x,y,alpha)
%Compute Bland-Altman plot characteristics.
%
%The method for calculating these values was taken from Bland&Altman
%(Lancet, 1986, i:307-310) http://dx.doi.org/10.1016/S0140-6736(86)90837-8
%
%The implementation of tinv is by Star Strider
% https://www.mathworks.com/matlabcentral/fileexchange/56500
%The implementation of alpha_to_Z is by Star Strider as well
% https://www.mathworks.com/matlabcentral/answers/45173#answer_55318

% http://web.archive.org/web/20190326204734/
% https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/
%  submissions/56500/versions/1/contents/tstat3.m  
%
% http://web.archive.org/web/20150317011457/
% http://www.mathworks.com/matlabcentral/answers/
%  45173-p-value-to-z-score#answer_55318

%Retrieve statistical functions.
persistent alpha_to_Z
if isempty(alpha_to_Z)
    if ifversion('>=',7.0,'Octave','>',3)
        %Anonymous functions were introduced in ML7.0, so we need to use eval to fool the syntax
        %checker in ML6.5.
        alpha_to_Z = eval('@(a) -sqrt(2) * erfcinv(a*2)'); % Get Z value given the alpha.
    else
        alpha_to_Z=inline('-sqrt(2) * erfcinv(a*2)','a'); %#ok<DINLN>
    end
end

%Compute main characteristics.
data=y(:)-x(:);data(isnan(data))=[];
n=numel(data);%Cases
mu=mean(data);%Mean
s=std(data);  %SD

%Compute statistical elements.
mu_ste=sqrt(s^2/n);     %Standard error of mean
loa_ste=sqrt(3*s^2/n);  %Standard error of limit of agreement
t=t_inv(alpha/2,n-1);   %t-value
z=alpha_to_Z(1-alpha/2);%Z value for two-tailed alpha

%Compute output parameters.
CI.mu=[mu-t*mu_ste mu+t*mu_ste];
loa=[mu-z*s mu+z*s];
CI.loa_lower=[loa(1)-(t*loa_ste) loa(1)+(t*loa_ste)];
CI.loa_upper=[loa(2)-(t*loa_ste) loa(2)+(t*loa_ste)];
end
function out=t_inv(a,v)
% This implementation of tinv is by Star Strider
%
% https://www.mathworks.com/matlabcentral/fileexchange/56500
%
% http://web.archive.org/web/20190326204734/
% https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/
%  submissions/56500/versions/1/contents/tstat3.m
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
persistent oct__tinv
if isOctave
    %Octave doesn't allow passing extra variables to fzero.
    if isempty(oct__tinv)
        %This uses eval to fool the syntax checker, but this requires some attention to make it
        %still work if the code is minified.
        oct__tdist2T=eval('@(t,v) (1-betainc(v/(v+t^2), v/2, 0.5));');
        oct__tdist1T=eval(['@(t,v) 1-(1-' var2str(oct__tdist2T) '(t,v))/2;']);
        oct__tinv=eval(['@(alpha,v) fzero(@(tval) (max(alpha,(1-alpha)) - ' ...
            var2str(oct__tdist1T) '(tval,v)), 5);']);
    end
    out=oct__tinv(a,v);
else
    out=fzero(@t_inv_reimplementation___local_fun,5,optimset('fzero'),a,v);
end
end
function out=t_inv_reimplementation___tdist2T(t,v)
out=1-betainc(v/(v+t^2), v/2, 0.5);
end
function out=t_inv_reimplementation___tdist1T(t,v)
out=1-(1-t_inv_reimplementation___tdist2T(t,v))/2;
end
function out=t_inv_reimplementation___local_fun(tval,a,v)
out=max(a,(1-a))-t_inv_reimplementation___tdist1T(tval,v);
end
function [success,options,ME]=BlandAltmanPlot_parse_inputs(var1,var2,varargin)
%Parse the inputs of the BlandAltmanPlot function.
% This function returns a success flag, the parsed options, and an ME struct. As input, the
% options should either be entered as a struct or as Name,Value pairs. Missing fields are filled
% from the default.

%Pre-assign outputs.
success=false;
options=struct;
ME=struct('identifier','','message','');

%Test the required inputs.
if numel(var1)==0 || numel(var1)~=max(size(var1)) || ~isnumeric(var1) ...
        || numel(size(var1))>2 %ndims must be <=2 for use in plot()
    ME.message='The first input is not a numeric vector.';
    ME.identifier='HJW:BlandAltmanPlot:incorrect_input';
    return
end
if numel(var2)==0 || numel(var2)~=max(size(var2)) || ~isnumeric(var2) ...
        || numel(size(var2))>2 %ndims must be <=2 for use in plot()
    ME.message='The second input is not a numeric vector.';
    ME.identifier='HJW:BlandAltmanPlot:incorrect_input';
    return
end
if numel(var1)~=numel(var2)
    ME.message='The first and second inputs do no match.';
    ME.identifier='HJW:BlandAltmanPlot:incorrect_input';
    return
end

%Set the defaults.
default.StoreToAppdata=true;  %Store output in the axis with setappdata.
default.alpha=0.05;           %Set default alpha value to 5%.
default.plot_x_mean=true;     %Use the mean as x coordinate.
default.plotCI=true;          %Add the CIs of the mean and LoAs to the plot.
default.GroupIndex=[];        %Don't split the plot.
default.AddDetailsText=true;  %Add the details to the plot as text.
default.TextDigitsDisplayed=4;%Number of digits or num2str FormatSpec.
default.xxyy=[];              %Axis range, leave empty or NaN to auto-determine.
default.TargetAxes=[];        %Axes where the plots and texts should be created.
default.Parent=[];            %Alternative option name for TargetAxes.

%The required inputs are checked, so now we need to return the default options if there are no
%further inputs.
if nargin==2
    options=default;
    success=true;
    options=BlandAltmanPlot_parse_inputs_ensureParentAxesExists(options);
    return
end

%Test the optional inputs.
struct_input=nargin==3 && isa(varargin{1},'struct');
NameValue_input=mod(nargin,2)==0 && ...
    all(cellfun('isclass',varargin(1:2:end),'char'));
if ~( struct_input || NameValue_input )
    ME.message=['The third input (options) is expected to be either a ',...
        'struct,',char(10),'or consist of Name,Value pairs.']; %#ok<CHARTEN>
    ME.identifier='HJW:BlandAltmanPlot:incorrect_input_options';
    return
end
if NameValue_input
    %convert the Name,Value to a struct
    for n=1:2:numel(varargin)
        try
            options.(varargin{n})=varargin{n+1};
        catch
            ME.message='Parsing of Name,Value pairs failed.';
            ME.identifier='HJW:BlandAltmanPlot:incorrect_input_NameValue';
            return
        end
    end
else
    options=varargin{1};
end
fn=fieldnames(options);
for k=1:numel(fn)
    curr_option=fn{k};
    item=options.(curr_option);
    ME.identifier=...
        ['HJW:BlandAltmanPlot:incorrect_input_opt_' lower(curr_option)];
    switch curr_option
        case 'alpha'
            if ~isnumeric(item) || numel(item)~=1 || ...
                    item<0 || item>1 || isnan(item)
                %Because of the short circuit or() the NaN test doesn't need to consider the case
                %of a vector, as the numel test will already become true, in which case it is not
                %evaluated.
                ME.message='The value of options.alpha is not a scalar between 0 and 1.';
                return
            end
        case {'StoreToAppdata','plotCI','AddDetailsText','plot_x_mean'}
            [passed,item]=test_if_scalar_logical(item);
            options.(curr_option)=item;
            if ~passed
                ME.message=['The value of options.',curr_option,' is not a logical scalar.'];
                return
            end
        case 'TextDigitsDisplayed'
            temp=false;
            %Test an integer input and a char input separately to keep the logical statements and
            %program flow more clear.
            if isnumeric(item)
                if numel(item)~=1 || double(item)~=round(double(item))
                    %The integer test could cause issues due to float rounding, but it should
                    %actually not be an issue.
                    %The integer test in its current form also catches NaN.
                    temp=true;
                end
            elseif isa(item,'char')
                try
                    %This will not catch all invalid FormatSpec inputs, but at least any errors are
                    %caught. If the user chooses to supply '%%d' as a FormatSpec, the function will
                    %run as expected (i.e. return '%d' on Octave and all Matlab releases except
                    %v6.5 where this returns an error).
                    temp_out=num2str(123.45,item);
                    if ~isa(temp_out,'char') || numel(temp_out)==0
                        temp=true;
                    end
                catch
                    temp=true;
                end
            else
                temp=true;
            end
            if temp
                ME.message='The value of options.TextDigitsDisplayed is not an integer scalar.';
                return
            end
        case 'xxyy'
            %NaN values will be filled with the data extent.
            if ~isnumeric(item) || (numel(item)~=4 && ~isempty(item)) || ...
                    ( ~isempty(item) && ...
                    ( item(2)<=item(1) || item(4)<=item(3) ) )
                ME.message=['The value of options.xxyy is neither ',...
                    'a valid 2D axis range, nor empty.'];
                return
            end
        case {'TargetAxes','Parent'}
            %It is difficult to test if a supplied variable is an axes object, as in old Matlab
            %(and in Octave) a handle is simply a double. A strategy with get (with try-catch)
            %should work.
            try
                if ~strcmp(get(item,'Type'),'axes') && ~isempty(item)
                    ME.message=['The value of options.',curr_option,...
                        ' is not an axes object or is not empty.'];
                    return
                end
            catch
                if ~isempty(item)
                    ME.message=['The value of options.',curr_option,...
                        ' is not an axes object or is not empty.'];
                    return
                end
            end
        case 'GroupIndex'
            if isempty(item),continue,end %If explicitly set to empty there is no need for checks.
            try
                if numel(item)~=numel(var1)
                    error('trigger')
                end
                x=zeros(1,max(item(:)));
                x(item)=1; %An indexing operation should work for valid input.
                if sum(x)==0 %At least 1 element should be marked.
                    error('trigger')
                end
            catch
                ME.message=['The values in GroupIndex must be valid array indices',...
                    char(10),'and the number of elements must match the data.']; %#ok<CHARTEN>
            end
        otherwise
            ME.message=sprintf('Name,Value pair not recognized: %s',curr_option);
            ME.identifier='HJW:BlandAltmanPlot:incorrect_input_NameValue';
            return
    end
end

%Fill any missing fields.
fn=fieldnames(default);
for k=1:numel(fn)
    if ~isfield(options,fn(k))
        options.(fn{k})=default.(fn{k});
    end
end

if isempty(options.TargetAxes)
    options=BlandAltmanPlot_parse_inputs_ensureParentAxesExists(options);
end
success=true;ME=[];
end
function options=BlandAltmanPlot_parse_inputs_ensureParentAxesExists(options)
%gca might result in unexpected behavior if called unintentionally
%
%When this function is called, the validity of the inputs is already checked, or the default was
%loaded (meaning both inputs are empty). The precedence is described below.
%
%behavior: (P='Parent' parameter, T='TargetAxes' parameter)
% isempty(P) &&  isempty(T) --> out=gca
%~isempty(P) &&  isempty(T) --> out=P
% isempty(P) && ~isempty(T) --> out=T
%~isempty(P) && ~isempty(T) --> out=P
P=options.Parent;T=options.TargetAxes;
ax_matrix={[],T;P,P};
ax=ax_matrix{2-isempty(P),2-isempty(T)};
if isempty(ax)
    ax=gca;
end
options.Parent=ax;
options.TargetAxes=ax;
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
% tf=ifversion(test,Rxxxxab)
% tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Output:
% tf       - If the current version satisfies the test this returns true.
%            This works similar to verLessThan.
%
% Inputs:
% Rxxxxab - Char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the
%           numeric version.
% test    - Char array containing a logical test. The interpretation of this is equivalent to
%           eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.6                                                         |%
%|  Date:    2021-03-11                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
        else
            L=ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf=NaN;return
            else
                v=v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%The char and string test are not case sensitive.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
%- 'enable' or 'disable'
%- 'enabled' or 'disabled'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    try
        states(end+1,:)=eval('{"on","off"}');
    catch
    end
end
isLogical=true;
try
    if isa(val,'char') || isa(val,'string')
        try val=lower(val);catch,end
    end
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end
function varargout=var2str(varargin)
%Analogous to func2str, return the variable names as char arrays, as detected by inputname.
%This returns an error for invalid inputs and if nargin~=max(1,nargout).
%
%You can use comma separated lists to create a cell array:
% out=cell(1,2);
% foo=1;bar=2;
% [out{:}]=var2str(foo,bar);
err_flag= nargin~=max(1,nargout) ;
if ~err_flag
    varargout=cell(nargin,1);
    for n=1:nargin
        try varargout{n}=inputname(n);catch,varargout{n}='';end
        if isempty(varargout{n}),err_flag=true;break,end
    end
end
if err_flag
    error('Invalid input and/or output.')
end
end
