classdef histogram2Polar < matlab.mixin.SetGet
    %histogram2Polar Bivariate histogram plot from polar coordinates / position vectors.
    % Esentially, the coordinates are converted to cartesian coordinates and fed to
    % HISTOGRAM2. A polar-coordinate grid is then drawn in order to mimic a polar plot.
    % The properties and methods of this class are intended to mimic a polar plot like it
    % is created by POLARAXES and properties with the same name, such as 'ThetaDir', mimic
    % the behavior of a real POLARAXES. Note that this has limitations and might not work
    % flawlessly. Properties may be specified via Name-Value pairs during plot
    % construction or set post-hoc using dot-notation or SET syntax. HISTOGRAM2 is called
    % with the parameters {'DisplayStyle', 'tile', 'ShowEmptyBins', 'off', 'EdgeColor',
    % 'none'}, which can be overridden by specifying the respective Name-Value pairs
    % during plot creation or by accessing the returned object's property 'Histogram'.
    % Optional Name-Value pairs of this function support those specified below and all
    % Name-Value pairs as accepted by HISTOGRAM2.
    %
    %
    % Syntax:
    %
    %   histogram2Polar(theta, rho);
    %   histogram2Polar(theta, rho, BinWidth);
    %   histogram2Polar(theta, rho, BinWidth, Name, Value);
    %   histogram2Polar(theta, rho, Name, Value);
    %   histogram2Polar(target, ___);
    %   h = histogram2Polar(___);
    %
    %
    % Input:
    %
    %   theta       Angular coordinate in radians. THETA and RHO must have the same size
    %               or any of them can be scalar.
    %
    %   rho         Radial coordinate. THETA and RHO must have the same size or any of
    %               them can be scalar.
    %
    %   BinWidth    Edge length of the square-sized bins in radial units. May be omitted
    %               or specified empty (default == 0.05). Must be specified if Name-Value
    %               pairs are given.
    %
    %   target      Handle to cartesian axes. Alternatively, the target axes can be
    %               specified as the 'Parent' parameter. If this syntax is used and the
    %               'Parent' parameter specified, the axes handle provided via TARGET is
    %               preferred over the one given as parameter value.
    %
    %
    % Name-Value pairs can be any pairs that are accepted by HISTOGRAM2 and additionally
    % the following (note that some Name-Value pairs from HISTOGRAM2 may be unreasonable):
    %
    %   'RTicks' -- Tick values of the radius axis
    %       -1 (default) | N-element vector of increasing values or []
    %       Specification as -1 leads to three equally spaced tick marks between the data
    %       limits. If empty, no ticks are drawn.
    %
    %   'RLim' -- Radius-axis limits
    %   	[0 1] (default) | Two-element numeric vector of increasing values
    %       Analogous to the POLARAXES property with the same name. If the radius data lie
    %       outside the initially specified limits, the limits are adjusted. Note,
    %       however, that changing the lower limit only affects the ticks (see note
    %       below).
    %
    %   'RAxisLocation' -- Location of radius-axis tick labels, in degrees
    %   	80 (default) | Scalar positive value
    %       Analogous to the POLARAXES property with the same name.
    %
    %   'ThetaTicks' -- Tick values of the angle axis, in degrees
    %   	-1 (default) | N-element vector of increasing values or []
    %       Specification as -1 leads to three equally spaced tick marks between the
    %       angle-axis limits. If empty, no ticks are drawn.
    %
    %   'ThetaLim' -- Angle-axis limits, in degrees
    %   	[0 360] | Two-element numeric vector of increasing values
    %       Analogous to the POLARAXES property with the same name.
    %
    %   'ThetaZeroLocation' -- Location of the zero reference axis
    %   	'top' (default) | 'left' | 'bottom' | 'right'
    %       Analogous to the POLARAXES property with the same name.
    %
    %   'ThetaDir' -- Direction of increasing angles
    %   	'counterclockwise' (default) | 'clockwise'
    %       Analogous to the POLARAXES property with the same name.
    %
    %   'ThetaAxisUnits' -- Units for angle values
    %       'degrees' (default) | 'radians'
    %       Analogous to the POLARAXES property with the same name.
    %
    %   'DisplayColorbar' -- Display color bar
    %   	true (default) | false
    %       Equals calling 'colorbar(_)', but has the advantage that the colorbar label
    %       displays the normalization method used by HISTOGRAM2 and that the handle to
    %       the colorbar is stored as a property.
    %
    %   'GridColor' -- Color of polar coordinate grid
    %   	[0.5 0.5 0.5] | valid MATLAB color specifier
    %
    %   'DisplayGridOnTop' -- Display polar-coordinate grid on top of histogram
    %   	false (default) | true
    %
    %   'Parent' -- Handle to the axes where the histogram is plotted
    %       current axes (default) | valid scalar handle to cartesian axes
    %       Alternatively, the target axes can be specified as a first input argument; in
    %       this case, specification of the 'Parent' parameter is ignored.
    %
    % All property values can be adjusted after plotting.
    %
    %
    % Output:
    % 	h       Handle to the created HISTOGRAM2POLAR object. Access its properties
    %           using dot notation or SET syntax.
    %
    %
    % Notable properties beside the input parameters (see DOC page for details):
    %
    % 	'Data', 'Histogram', 'Colorbar', 'GridEdge', 'RAxisGrid', 'RAxisLabels',
    % 	'ThetaAxisGrid', 'ThetaAxisLabels'
    %
    %   
    % Methods:
    %
    %   fewerbins, morebins
    %
    %
    %
    % Example:
    % theta = randn(3000, 1)/5;
    % rho = normalize(randn(size(theta)), 'range');
    % figure, polarplot(theta, rho, '.');
    % figure, histogram2Polar(theta, rho);
    % figure, histogram2Polar(theta, rho, 0.1, 'ThetaTicks', 0:45:360, 'thetaZeroLocation', 'right', 'edgeColor', 'k', 'Normalization', 'pdf');
    %
    %
    %
    % Notes:
    % - The radius-axis lower limit can be changed, but the origin will always be at
    %   (0,0). Implementing an actual limits-change would be rather difficult (I suppose)
    %   because the radius data, ticks and labels would need scaling. Maybe this will be a
    %   feature of a future version.
    % - Minor grid and tick values are currently unsupported.
    % - Displaying angle units in radians does not format the angles as multiples of pi,
    %   as opposed to the behavior of POLARAXES. This is unfortunate. 
    %       
    %
    %
    % ---Author: Frederick Zittrell
    %
    % See also histogram2 polaraxes colormap set
    
    properties
        % N-by-2 matrix where the first column contains angle data and the second column contains radius data. Translates to the same property of HISTOGRAM2 and adjusts the plot if changed.
        Data double
        
        % Positive scalar. Edge length of square-sized bins, in radius data. Adjustable after plotting.
        BinWidth(1,1) double {mustBeNonnegative} = histogram2Polar.defBinWidth
        
        % Handle to the AXES where the histogram is plotted.
        Parent(1,1) matlab.graphics.axis.Axes {histogram2Polar.mustBeValidAx} = gca
        
        % Angle-axis tick values (in degrees), same behavior as for POLARAXES.
        ThetaTicks double = histogram2Polar.defThetaTicks
        
        % 'degrees' | 'radians'. Angle axis units, same behavior as for POLARAXES.
        ThetaAxisUnits char {histogram2Polar.mustBeValidUnit} = histogram2Polar.defThetaUnits
        
        % Angle-axis limits (in degrees), same behavior as for POLARAXES.
        ThetaLim double {histogram2Polar.mustBeValidLim} = histogram2Polar.defThetaLim
        
        % 'top' | 'left' | 'right' | 'bottom'. Location of the zero reference axis, same behavior as for POLARAXES.
        ThetaZeroLocation char {histogram2Polar.mustBeValidLoc} = histogram2Polar.defThetaZeroLoc
        
        % 'counterclockwise' | 'clockwise'. Direction of increasing angles, same behavior as for POLARAXES.
        ThetaDir char {histogram2Polar.mustBeValidThetaDir} = histogram2Polar.defThetaDir
        
        % Radius-axis tick values, same behavior as for POLARAXES.
        RTicks double = histogram2Polar.defRTicks
        
        % Radius-axis limits, same behavior as for POLARAXES.
        RLim double {histogram2Polar.mustBeValidLim} = histogram2Polar.defRLim
        
        % Location of radius-axis tick labels (in degrees), same behavior as for POLARAXES.
        RAxisLocation (1,1) double {mustBeNumeric} = histogram2Polar.defRAxisLoc
        
        % Color of polar-coordinate grid lines.
        GridColor = histogram2Polar.defGridColor
        
        % Boolean, controls whether the polar-coordinate grid is displayed on top of the histogram or beneath.
        DisplayGridOnTop double {mustBeNumericOrLogical} = histogram2Polar.defDisplayGridOnTop
        
        % Boolean, controls whether a color bar is displayed.
        DisplayColorbar double {mustBeNumericOrLogical} = histogram2Polar.defDisplayColorbar
    end
    
    
    properties (SetAccess = private)
        % handles to graphic objects
        
        % Histogram2 object. Access histogram properties via this handle. See also <a href="matlab: doc histogram2">histogram2</a>
        Histogram
        
        % Color bar handle. Use this to change color bar properties; only valid if 'DisplayColorbar' is TRUE.
        Colorbar
        
        % Handle to the grid circle at the plot edge. Shares 'Color' and 'LineWidth' with the other grid lines.
        GridEdge
        
        % Handle to the radius-axis tick-circles. Shares 'Color' and 'LineWidth' with the other grid lines.
        RAxisGrid
        
        % Handle to the radius-axis labels, which are separate TEXT objects and share font properties.
        RAxisLabels
        
        % Handle to the angle-axis grid lines. Shares 'Color' and 'LineWidth' with the other grid lines.
        ThetaAxisGrid
        
        % Handle to the angle-axis labels, which are separate TEXT objects and share font properties.
        ThetaAxisLabels
    end
    
    
    properties (Dependent, Access = private)
        thDirFac    % Stores the angle multiplication sign for correct data and label placement depending on 'ThetaDir'.
    end
    
    properties (Constant, Access = private)
        % default property values and validation functions
        
        defBinWidth = 0.05
        
        defThetaTicks = -1
        defRTicks = -1
        
        defThetaUnits = 'degrees'
        mustBeValidUnit = @(c) mustBeMember(lower(c), {'degrees', 'radians'})
        
        defThetaDir = 'counterclockwise'
        mustBeValidThetaDir = @(c) mustBeMember(lower(c), {'counterclockwise', 'clockwise'})
        
        defThetaZeroLoc = 'top'
        mustBeValidLoc = @(c) mustBeMember(lower(c), {'top', 'left', 'bottom', 'right'})
        
        defThetaLim = [0 360]
        defRLim = [0 1]
        mustBeValidLim = @(v) validateattributes(v, {'numeric'}, {'increasing', 'numel', 2, 'vector'})
        
        mustBeValidAx = @(h) assert(isscalar(h) && isa(h, 'matlab.graphics.axis.Axes') && isvalid(h), ...
            'Target axes must be a valid axes handle.')
        
        defRAxisLoc = 80
        
        defGridColor = [1 1 1] * 0.5
        
        defDisplayColorbar = true
        
        defDisplayGridOnTop = false
        
        %
        circleResolution = 0.02 % drawing resolution of RTick circles (in radians)
    end
    
    
    methods
        function self = histogram2Polar(theta, rho, BinWidth, varargin)
            %histogram2Polar Constructor: histogram2Polar(theta, rho, binWidth, Name, Value)
            %
            %   obj = histogram2Polar(theta, rho, binWidth, Name, Value)
            
            %% input parsing           
            assert(nargin >= 2, 'Not enough input arguments');
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            
            ip.addRequired('theta', @mustBeNumeric);
            ip.addRequired('rho', @mustBeNumeric);
            
            ip.addOptional('BinWidth', self.defBinWidth);
            
            % defining default values here _and_ in the property definition is obviously
            % redundant. Note that the property-value validation is only defined in the
            % property definition, while the default values are specified in the defintion
            % and here in the parser. I like defining the validations in the defintions
            % because the structure is much cleaner than doing this in the parser and
            % because the generated error messages when a validation function fails are
            % very nice. However, to use the parser the way it is intended (automatically
            % setting specified and default values), the default values need to be given
            % to the parser as well. Inputs THETA and RHO are different because these are
            % not assigned property values.
            ip.addParameter('RTicks', self.defRTicks);
            ip.addParameter('RLim', self.defRLim);
            ip.addParameter('RAxisLocation', self.defRAxisLoc);
            ip.addParameter('ThetaTicks', self.defThetaTicks);
            ip.addParameter('ThetaLim', self.defThetaLim);
            ip.addParameter('ThetaZeroLocation', self.defThetaZeroLoc);
            ip.addParameter('ThetaDir', self.defThetaDir);
            ip.addParameter('ThetaAxisUnits', self.defThetaUnits);
            ip.addParameter('GridColor', self.defGridColor);
            ip.addParameter('DisplayColorbar', self.defDisplayColorbar);
            ip.addParameter('DisplayGridOnTop', self.defDisplayGridOnTop);
            ip.addParameter('Parent', gca);
            
            if nargin < 3 % easy
                ip.parse(theta, rho);
            else % hard: facilitating optional TARGET specification as first argument
                if isscalar(theta) && isa(theta, 'matlab.graphics.axis.Axes') && isvalid(theta)
                    % first input is axes -> shift input variable contents to the left and
                    % use first input as 'Parent' value; note that by appending it to
                    % VARARGIN, this causes the first argument to have preference over the
                    % given parameter value if 'Parent' is also specified via Name-Value
                    % pair input
                    varargin = [varargin, {'Parent', theta}];
                    theta = rho;
                    rho = BinWidth;
                    
                    % now, either there were more arguments specified, then BINWIDTH is
                    % the first element of VARARGIN -> shift contents
                    if nargin > 3
                        BinWidth = varargin{1};
                        varargin = varargin(2:end);
                        
                    else
                        % or no additional arguments were given -> manually assign default
                        % value to BINWIDTH
                        BinWidth = self.defBinWidth;
                    end
                end
                
                % BINWIDTH may also be specified as [] -> assign default value
                if isempty(BinWidth), BinWidth = self.defBinWidth; end
                
                % finally, parse everything; notes regarding BINWIDTH: If it is omitted,
                % it is assigned the default value by the parser because it is not passed
                % to it (see above PARSE call). If it is omitted but Name-Value pairs are
                % given, the below PARSE call still works because then, BINWIDTH is a
                % parameter name, being recognized by the parser as such, and the parser
                % returns the default value for BINWIDTH.
                ip.parse(theta, rho, BinWidth, varargin{:});
            end
            
            theta = ip.Results.theta(:);
            rho = ip.Results.rho(:);
            assert(numel(rho) == numel(theta), ...
                'Inputs THETA and RHO must have equal numbers of elements.');
            self.Data = [theta, rho];
            
            % assign parameters to properties; exclude THETA and RHO from parser
            % parameter-list
            params = setdiff(ip.Parameters, {'theta', 'rho'});
            for iParam = 1:numel(params)
                cParam = params{iParam};
                self.(cParam) = ip.Results.(cParam);
            end
            
            % pass unmatched input to HISTOGRAM2: Create cell matrix with input that could
            % not be matched by the parser, then transform the matrix to reconstruct the
            % name-value pairs
            hist2Args = [fieldnames(ip.Unmatched)'; struct2cell(ip.Unmatched)'];
            hist2Args = hist2Args(:)';
            
            
            % default dynamic values
            if self.RLim(2) < max(rho)
                self.RLim = [0 ceil(max(rho))]; end
            if isscalar(self.RTicks) && self.RTicks < 0
                self.RTicks = linspace(self.RLim(1), self.RLim(2), 5);
            end
            
            if isscalar(self.ThetaTicks) && self.ThetaTicks == -1
                self.ThetaTicks = linspace(self.ThetaLim(1), self.ThetaLim(2), 5);
            end
            
            axH = self.Parent;
            
            %% 
            hold(axH, 'on');
            
            % polar grid is drawn implicitly during tick-value setting
            self.drawEdge;
            
            %% create histogram
            % coordinate transformation
            [x, y] = pol2cart(theta * self.thDirFac, rho);
            
            edges = self.getEdges;
            
            try
                self.Histogram = histogram2(axH, x, y, edges, edges, ...
                    'DisplayStyle', 'tile', 'ShowEmptyBins', 'off', 'edgecolor', 'none', ...
                    hist2Args{:});
            catch ME
                if strcmp(ME.identifier, 'MATLAB:unrecognizedStringChoice')
                    % handle exception of non-matching parameter name
                    
                    validParams = ip.Parameters; % compose parameter names
                    validParams = setdiff(validParams, {'theta', 'rho'});
                    validParams = cellfun(@(c) ['''',c,''''], validParams, 'UniformOutput', false);
                    exMsg = sprintf('Unrecognized Name-Value pair input\nValid parameters for %s:\n\n%s\n', ...
                        '<a href="matlab: matlab.internal.language.introspective.errorDocCallback(''histogram2Polar'')">histogram2Polar</a>', ...
                        strjoin(validParams, ', '));
                    
                    % create exception containing the valid parameters of histogram2Polar,
                    % add the caught exception as cause and throw it. This exception first
                    % lists the valid properties and then outputs the histogram2 exception
                    % message. Neat.
                    ex = MException('HISTOGRAM2POLAR:unrecognizedParameterName', exMsg);
                    ex = ex.addCause(ME);
                    ex.throw;
                else, ME.rethrow;
                end
            end
           
            % if histogram is deleted, the respective instance of this class is deleted as
            % well
            self.Histogram.DeleteFcn = @(~,~) self.delete;
            
            %%
            colormap(self.Parent, jet);
            self.drawColorbar;
            
            
            % axis operations
            axis(axH, 'off');
            axis(axH, 'equal');
            axis(axH, 'tight');
            xlim(axH, [-self.RLim(2), self.RLim(2)]);
            ylim(axH, [-self.RLim(2), self.RLim(2)]);
            
            % link axes properties to the same properties of the polar grid, so changing
            % the axes properties updates the polar grid and changing the same property of
            % one grid element adjusts the property of the other elements
            axH.UserData.propLinkLabelsTheta = linkprop( ...
                [axH.YAxis; self.ThetaAxisLabels], {'FontSize', 'FontWeight', 'FontName'});
            axH.UserData.propLinkLabelsRho = linkprop( ...
                [axH.XAxis; self.RAxisLabels], {'FontSize', 'FontWeight', 'FontName'});
            
            axH.UserData.propLinkLines = linkprop( ...
                [axH.XAxis; axH.YAxis; self.RAxisGrid; self.ThetaAxisGrid; self.GridEdge], ...
                {'Color', 'LineWidth'});
            axH.XAxis.Color = self.GridColor;
            axH.YAxis.Color = self.GridColor;
            
            % link parent axes of all graphic objects -- this is probably fully
            % unnecessary.
            axH.UserData.propLinkParent = linkprop([self.Histogram; ...
                self.RAxisGrid; self.GridEdge; self.RAxisLabels; ...
                self.ThetaAxisGrid; self.ThetaAxisLabels], 'Parent');
        end % constructor
        
        %%
        function set.ThetaTicks(self, tt)
            %ThetaTicks
            self.ThetaTicks = tt;
            if ~(isscalar(self.ThetaTicks) && self.ThetaTicks == -1)
                self.drawThetaTicks;
                self.drawThetaLabels;
            end
        end
        
        
        function set.ThetaLim(self, lim)
            %ThetaLim Only affects the grid and ticks. To exclude data points, change
            %'Data' directly.
            
            self.ThetaLim = lim;
            
            % -1 case occurs during construction with default THETATICK value. I'm
            % terribly sorry if you really want to have a single tick at -1 -- consider
            % specifying the tick value as 359.
            if ~isempty(self.ThetaTicks) && ~(isscalar(self.ThetaTicks) && self.ThetaTicks ~= -1)
                % adjust ticks: remove ticks outside limits
                oldTicks = mod(self.ThetaTicks, 360);
                ticks = oldTicks(oldTicks >= lim(1) & oldTicks <= lim(2));
                
                self.ThetaTicks = ticks;
                self.drawEdge;
                self.drawRTicks;
            end
        end
        
        
        function set.RTicks(self, rt)
            %RTicks
            self.RTicks = rt;
            if ~(isscalar(self.RTicks) && self.RTicks < 0)
                self.drawRTicks;
                self.drawRLabels;
            end
        end
        
        
        function set.BinWidth(self, bw)
            %BinWidth Adjusts the bin edge length.
            if bw ~= self.BinWidth
                self.BinWidth = bw;
                
                if ~isempty(self.Histogram) && isvalid(self.Histogram)
                    edges = self.getEdges;
                    set(self.Histogram, 'XBinEdges', edges, 'YBinEdges', edges);
                end
            end
        end
        
        
        function set.Data(self, d)
            %Data   Changes the corresponding histogram data. The new DATA must be a
            %N-by-2 matrix where the first column contains the angles and the second
            %column contains the radius, analogous to the behavior of HISTOGRAM2.
            self.Data = d;
            
            if ~isempty(self.Histogram) && isvalid(self.Histogram)
                [x, y] = pol2cart(d(:,1), d(:,2));
                self.Histogram.Data = [x, y];
            end
        end
        
        
        function set.DisplayGridOnTop(self, tf)
            %DisplayGridOnTop
            if tf ~= self.DisplayGridOnTop % only operate on change
                self.DisplayGridOnTop = tf;
                
                % get grid line handles
                lineH = vertcat(self.ThetaAxisGrid, self.RAxisGrid, self.GridEdge);
                
                if ~isempty(lineH)
                    lineH = lineH(isvalid(lineH));
                    
                    if tf, z = 0;
                    else,  z = -0.1; end
                    
                    % adjust Z-coordinate
                    arrayfun(@(h) set(h, 'ZData', repmat(z, size(h.XData))), lineH);
                end
            end
        end
        
        
        function set.Parent(self, p)
            % Sets the same Parent to the histogram. Not sure whether this is even
            % beneficial.
            self.Parent = p;
            if ~isempty(self.Histogram) && isvalid(self.Histogram)
                self.Histogram.Parent = p; end
        end
        
        
        function set.DisplayColorbar(self, tf)
            %DisplayColorbar
            self.DisplayColorbar = tf;
            if tf
                self.drawColorbar;
            elseif ~isempty(self.Colorbar) && isvalid(self.Colorbar)
                delete(self.Colorbar);
            end
        end
        
        
        function set.ThetaZeroLocation(self, loc)
            %ThetaZeroLocation Mimics the behavior of the same property of a POLARAXES.
            self.ThetaZeroLocation = loc;
            switch lower(loc)
                case 'top',    viewAz = -90;
                case 'left',   viewAz = -180;
                case 'bottom', viewAz = 90;
                case 'right',  viewAz = 0;
            end
            if isvalid(self.Parent)
                view(self.Parent, viewAz, 90); end
        end
        
        
        function set.ThetaAxisUnits(self, unit)
            %ThetaAxisUnits Mimics the behavior of the same property of a POLARAXES.
            if ~strcmpi(unit, self.ThetaAxisUnits)
                self.ThetaAxisUnits = unit;
                self.drawThetaLabels;
            end
        end
        
        
        function set.ThetaDir(self, tdir)
            %ThetaDir Mimics the behavior of the same property of a POLARAXES.
            toggle = ~strcmpi(self.ThetaDir, tdir);
            self.ThetaDir = tdir;
            
            if toggle
                % mirror data if THETADIR is _changed_
                self.Data(:,1) = self.Data(:,1) * -1;
                
                % update grid
                self.drawThetaTicks;
                self.drawThetaLabels;
                self.drawRTicks;
                self.drawRLabels;
            end
        end
        
                
        function thDirFac = get.thDirFac(self) % private
            %thDirFac Returns 1 or -1 based on 'ThetaDir'; used as a multiplication factor
            %for angle calculations to match 'ThetaDir'. For reverse theta direction,
            %angles are multiplied with -1; this also affects labels.
            if strcmpi(self.ThetaDir, 'counterclockwise'), thDirFac =  1;
            else,                                          thDirFac = -1; end
        end
        
        
        % morebins/fewerbins
        function morebins(self)
            %morebins Increases the number of bins.
            %
            %   obj.morebins;
            %
            % See also histogram2/morebins
            if ~isempty(self.Histogram)
                self.Histogram.morebins; end
        end
        function fewerbins(self)
            %fewerbins Decreases the number of bins.
            %
            %   obj.fewerbins;
            %
            % See also histogram2/fewerbins
            if ~isempty(self.Histogram)
                self.Histogram.fewerbins; end
        end
    end % public
    %%
    methods (Access = private)
        function edges = getEdges(self)
            %getEdges Compose histogram edges; this procedure assures that the bins are
            %symmetricrelative to the plot center, i. e., the innermost four bins meet at
            %(0,0).
            
            edgesPos = 0 : self.BinWidth : self.RLim(2);
            edgesNeg = -1 * edgesPos(end:-1:2);
            edges = [edgesNeg, edgesPos];
        end
        
        function ticksRad = uniqueThetaTicksRad(self)
            %uniqueThetaTicksRad Returns the unique tick marks in radians. Useful because
            %'ThetaLim' may contain the same circular value.
            ticksRad = deg2rad(sort(unique(mod(self.ThetaTicks, 360)), 'ascend'));
        end
        
        
        function drawThetaTicks(self)
            %drawThetaTicks
            if ~isempty(self.ThetaAxisGrid) && all(isvalid(self.ThetaAxisGrid))
                delete(self.ThetaAxisGrid); end
            
            theta = self.uniqueThetaTicksRad;
            
            % compile tick coordinate data: 2-by-N matrix where N is the number of ticks.
            % The upper row is zeros for both angle and radius because tick lines start at
            % the origin, the lower row is the upper RLim for the radius and the
            % respective tick value for the angle.
            thetaTicksRho = [zeros(size(theta)); repmat(self.RLim(2), size(theta))];
            theta = [zeros(size(theta)); theta];
            self.ThetaAxisGrid = self.drawGrid(theta, thetaTicksRho);
        end
        
        
        function drawThetaLabels(self)
            %drawThetaLabels
            if ~isempty(self.ThetaAxisLabels) && all(isvalid(self.ThetaAxisLabels))
                delete(self.ThetaAxisLabels); end
            
            theta = self.uniqueThetaTicksRad;
            
            % add offset to radial coordinate so labels are outside the plot
            thetaLabelRho = repmat(self.RLim(2), size(theta)) + range(self.RLim) * 0.15;
            [thetaLabelX, thetaLabelY] = pol2cart(theta', thetaLabelRho');
            
            % label text
            if strcmpi(self.ThetaAxisUnits, 'degrees')
                  num2Lbl = @(v) [num2str(rad2deg(v)),'°'];
            else, num2Lbl = @num2str;
            end
            labels = arrayfun(num2Lbl, mod(theta * self.thDirFac, 2*pi), ...
                'UniformOutput', false);
            
            self.ThetaAxisLabels = text(self.Parent, thetaLabelX, thetaLabelY, labels, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        
        function drawRTicks(self)
            %drawRTicks
            if ~isempty(self.RAxisGrid)
                delete(self.RAxisGrid); end
            
            if ~isempty(self.RTicks)
                thetaLim = deg2rad(self.ThetaLim);
                circTheta = (thetaLim(1) : self.circleResolution : thetaLim(2))';
                rho = self.RTicks;
                
                self.RAxisGrid = self.drawGrid(repmat(circTheta, 1, numel(rho)), rho);
            end
        end
        
        
        function drawRLabels(self)
            %drawRLabels
            if ~isempty(self.RAxisLabels)
                delete(self.RAxisLabels); end
            
            rho = self.RTicks;
            theta = repmat(deg2rad(self.RAxisLocation * self.thDirFac), size(rho));
            [rhoLabelsX, rhoLabelsY] = pol2cart(theta, rho);
            
            % convert tick values to char
            labels = arrayfun(@num2str, rho, 'UniformOutput', false);
            
            self.RAxisLabels = text(self.Parent, rhoLabelsX, rhoLabelsY, labels, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        
        function drawEdge(self)
            %drawEdge Draws the outer plot edge.
            if ~isempty(self.GridEdge)
                delete(self.GridEdge); end
            
            thetaLim = deg2rad(self.ThetaLim);
            circTheta = (thetaLim(1) : self.circleResolution : thetaLim(2))';
            self.GridEdge = self.drawGrid(circTheta, self.RLim(2));
        end
        
        
        function drawColorbar(self)
            if ~isempty(self.Histogram) && isvalid(self.Histogram) && self.DisplayColorbar
                self.Colorbar = colorbar(self.Parent);
                self.Colorbar.Label.String = self.Histogram.Normalization;
            end
        end
        
        
        function h = drawGrid(self, thetaRad, rho)
            %drawGrid   Helper function that draws grid elements.
            
            % if the lines and the histogram have the same Z-coordinates, the lines will
            % be on top of the histogram. To change this, the lines' Z-coordinates are
            % shifted below zero. (No, using UISTACK does not work.)
            if self.DisplayGridOnTop, z = 0;
            else,                     z = -0.1; end
            
            [circX, circY] = pol2cart(thetaRad, rho);
            h = line(self.Parent, circX, circY, repmat(z, size(circX)), 'Color', self.GridColor);
        end
    end % private
end