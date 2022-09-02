function MPUI

%MPUI is a graphical user interface to explore the geometric properties of moiré
% patterns (MPs) in the context of 2-dimensional materials and atomic lattices.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLEASE CONSULT THE README FILE FOR INSTRUCTIONS AND FURTHER DETAILS!
%
% Authors:
% Maxime Le Ster, Tobias Maerkl, Simon Brown
% The MacDiarmid Institute for Advanced Materials and Nanotechnology &
% School of Physical and Chemical Sciences
% University of Canterbury, Private Bag 4800, Christchurch 8140, New Zealand
%
%
%
% REMARKS ABOUT CODE:
%	1) Representation of vectors
%	All vectors are expressed as complex numbers, with the Real/Imaginary part
%	denoting the x/y component; this makes it easy to calculate the differences,
%	as well as carrying out rotations (a simple complex phase factor).
%
%	2) Units used in this program
% - all lengths are measured in Ångström;
% - inverse lenghts (k vectors) in Å^-1
%	- angle variables are labelled with their units in the code; however the user
%		is exposed only to angles that are expressed in degree!!
% - the positions of the second atom(s) in the basis are expressed as relative
%		coordinates, i.e. they must be between 0 and 1
%	
%
% 3) Convention
% We are using the convention 2*pi = 1 (crystallographer's convention) when
% going back and forth between real and reciprocal space!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CHECKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCREEN RESOLUTION CHECK
g = get(0);
if g.ScreenPixelsPerInch ~= 96
	resBoxString = {strcat('Warning: Your screen pixel density is ', char(160), num2str(g.ScreenPixelsPerInch), char(160)),...
		'i.e. it deviates from 96 dpi.', '',...
		'As a consequence the UI might not be displayed as intended.',...
		'Try changing your OS settings if you experience problems!',...
		'Also, please consult the Readme file for known issues.'};
	resBoxPosition = [(g.ScreenSize(3)-400)/2, (g.ScreenSize(4)-180)/2, 400, 180];
	resBox = figure('Position', resBoxPosition);
	resBoxButton = uicontrol(resBox, 'Position', [150 10 100 30], 'String', 'OK', ...
              'Callback', 'uiresume(gcbf);close(gcbf)');
	resBoxMessage = uicontrol('Style', 'text', 'Position', [5 42 390 130], 'String', resBoxString);
	uiwait(gcf);
end

% VERSION CHECK - R2016A OR NEWER
if verLessThan('matlab','9.0')
	VersionBoxString = {...
		'This program was developed and tested on Matlab R2016a & R2017a.',...
		'Your version seems to be older than that.',...
		'The program may not work or behave not as intended. Sorry!'};
	VersionBoxPosition = [(g.ScreenSize(3)-400)/2, (g.ScreenSize(4)-180)/2, 400, 180];
	VersionBox = figure('Position', VersionBoxPosition);
	VersionBoxButton = uicontrol(VersionBox, 'Position', [150 10 100 30], 'String', 'OK', ...
              'Callback', 'uiresume(gcbf);close(gcbf)');
	VersionBoxMessage = uicontrol('Style', 'text', 'Position', [5 42 390 130], 'String', VersionBoxString);
	uiwait(gcf);
end

% STARTUP DIALOG
StartUpString = {'Thank you for using our moiré pattern UI!','',...
	'By using this program you agree to cite the following research article',...
	'if this program contributes to a publication which you are (co-)authoring:',...
	'',...
	'M. Le Ster, T. Märkl & S. Brown, "Moiré Patterns: A Simple Analytical Model"',...
	'2019, 2D Materials 7, 011005'};
StartUpBoxPosition = [(g.ScreenSize(3)-400)/2, (g.ScreenSize(4)-180)/2, 400, 180];
StartUpBox = figure('pos', StartUpBoxPosition);
StartUpBoxButton = uicontrol('Position', [150 10 100 30], 'String', 'OK', ...
              'Callback', 'uiresume(gcbf);close(gcbf)');
StartUpBoxMessage = uicontrol('Style', 'text', 'Position', [5 42 390 130], 'String', StartUpString);
uicontrol('Style', 'pushbutton', 'Position', [85 50 230 25],...
	'String', 'https://dx.doi.org/10.1088/2053-1583/ab5470',...
	'ForegroundColor','blue',...
	'Callback', @myweb);
uiwait(gcf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iNan = NaN*exp(1i); % convenient complex NaN

[geom, rot, rsp, ksp, MP, thetaPanel, table, vispn, geopn] = createDataStructures;

h = addUIcomponents;


% The following commands deactivate several pushbuttons of the Matlab toolbar.
% Using them could break the UI or at least severely mess up the layout. Comment
% the following lines of you want full access to the toolbar.
toolBarList = {'Show Plot Tools and Dock Figure', 'Link Plot', 'Hide Plot Tools', 'Insert Legend',...
	'Insert Colorbar', 'Rotate 3D', 'Edit Plot'};

for i=1:length(toolBarList)
	foundToolBarButton = findall(h.fig, 'ToolTipString', toolBarList{i});
	set(foundToolBarButton, 'Visible', 'Off');
end
%uipushtool buttons
% (1) Show Plot Tools
% (2) Hide Plot Tools
% (3) Print Figure
% (4) Save Figure
% (5) Open File
% (6) New Figure
% uitoggletool buttons
% (1) Insert Legend
% (2) Insert Colorbar
% (3) Data Cursor
% (4) Rotate 3D
% (5) Pan
% (6) Zoom Out
% (7) Zoom In
% (8) Edit Plot
% other
% 'Show Plot Tools and Dock Figure'
% 'Link Plot'


recallGeometry('bilayer graphene (default)'); % {'bilayer graphene (default)', 'G on h-BN', 'a-Bi on MoS2', 'a-Sb on a-Bi'}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define data structures that are shared among
%	most of the panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [geom, rot, rsp, ksp, MP, thetaPanel, table, vispn, geopn] = createDataStructures
% "geom": struct that holds the information about the GEOMETRY of the OL and
% UL - and their geometry ONLY! The twist angle theta is a separate variable in
% the struct "rot".
	geom.ULUnitCellR1 = [];
	geom.ULUnitCellR2 = [];
	geom.ULUnitCellAngleDeg = [];
	geom.ULUnitCellAngleRad = [];
	geom.OLUnitCellR1 = [];
	geom.OLUnitCellR2 = [];
	geom.OLUnitCellAngleDeg = [];
	geom.OLUnitCellAngleRad = [];
	geom.ULSecondAtomR1 = []; % 2nd atom, relative coordinates between 0 and 1
	geom.ULSecondAtomR2 = [];
	geom.OLSecondAtomR1 = [];
	geom.OLSecondAtomR2 = [];
	geom.untwistedOL = []; % full, but untwisted (only global angle) OL lattice
	geom.untwistedOLBaseV1 = [];
	geom.untwistedOLBaseV2 = [];
	geom.ULBaseV1 = []; % complex numbers that represent the actual base vectors
	geom.ULBaseV2 = [];
	geom.OLBaseV1 = [];
	geom.OLBaseV2 = [];

% "rot": struct to hold the two rotation angles (global and twist)
	rot.TwistAngleDeg = [];
	rot.TwistAngleRad = [];
	rot.GlobalAngleDeg = []; % typically zero unless trying to match experimental data
	rot.GlobalAngleRad = [];

% "rsp": everything needed to calculate and display the Real-Space Panel
% (from now on "rsp")
	rsp.ULnUnitCells = 10; % positive number to determine the number of UL unit cells to be displayed, defaults to 10
	rsp.OLnUnitCells = 10; % same for OL
	rsp.ULSecondAtomIsVisible = false;
	rsp.OLSecondAtomIsVisible = false;
	rsp.ULCoords = [];	% complex 2-dim NxN matrix to hold the final, rotated, atomic coordinates
	rsp.OLCoords = [];	% same
	rsp.defaultScale = []; % default size for RSP display, is based on unit cells
	rsp.currentScale = [];
	rsp.markersize_1 = 18;
	rsp.markersize_2 = 20;
	rsp.resetScaleFlag = true; % this flag indicates whether the real space view should be reset to a default zoom

% "ksp": everything needed to calculate and display the Reciprocal/K-Space Panel
% (from now on "ksp")
	% display properties of the panel
	ksp.defaultScale = [-0.5 0.5 -0.5 0.5]; % Å^-1; this is the default size for KSP display;
	ksp.currentScale = []; 
	ksp.resetScaleFlag = true;
	ksp.annotationFlag = false;
	ksp.markerSize1 = 20; % marker size for the difference vectors Kmnpq
	ksp.markerSize2 = 40; % - " -	: OL and UL dots
	% quantities used in calculation of coordinates etc:
	ksp.cutoff = [];
	ksp.CutoffLimits = [0.01, 5]; % these are dummy default values; the program overrides them as soon as it loads an example
	ksp.ULnUnitCells = 4; % Number of k-space unit cells, default values
	ksp.OLnUnitCells = 4;
	ksp.ULk1 = []; % "lattice constants", only the LENGTH of k vectors (real numbers)
	ksp.ULk2 = [];
	ksp.OLk1 = [];
	ksp.OLk2 = [];
	ksp.ULBaseK1 = []; % k-space base vectors (complex numbers)
	ksp.ULBaseK2 = [];
	ksp.OLBaseK1 = [];
	ksp.OLBaseK2 = [];
	ksp.ULCoords = []; % these 2 matrices contain the full reciprocal lattices,
	ksp.OLCoords = []; % their size is based on ksp.ULnUnitCells & ksp.OLnUnitCells
	ksp.U = []; % ksp.U and ksp.O are trimmed from "Coords" to those values
	ksp.O = []; % smaller than the cut-off and excluding (0,0)
	ksp.allDiffVectorsK = [];	% finally, this is the variable to contain all the difference vectors between the two reciprocal lattices "OL-UL".

% "MP": struct to contain the actual MP properties, ie period and angle
	MP.l1 = []; % lambda1,2,3
	MP.l2 = [];
	MP.l3 = [];
	MP.angleDeg1 = [];
	MP.angleRad1 = [];
	MP.angleDeg2 = [];
	MP.angleRad2 = [];
	MP.angleDeg3 = [];
	MP.angleRad3 = [];
	MP.kList = []; % list for the resulting 3 shortest k vectors; used by MP table
	
% "thetaPanel": struct for lambda/delta vs (theta) panel
	thetaPanel.m = 1;
	thetaPanel.n = 0;
	thetaPanel.p = 1;
	thetaPanel.q = 0;
	thetaPanel.showflag = true; % This is needed for the startup phase, when the full geometry info is not yet there, to prevent the panel from being updated
	
% some UI element structs: "geopn" (geometry panel), "vispn" (visibility options),
% and "table" (output table for 3 sets of moire fringes)
	table = [];
	vispn = [];
	geopn = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lots of 'set' (and a few 'get') functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGeometryUL(R1, R2, w)
	geom.ULUnitCellR1 = R1;
	geom.ULUnitCellR2 = R2;
	geom.ULUnitCellAngleDeg = w;
	geom.ULUnitCellAngleRad = w*pi/180;
	
	% make calls to necessary update functions
	updateGeometryFields;
	if thetaPanel.showflag
		updateThetaPanel;
	end
end

function setGeometryOL(R1_, R2_, w_)
	geom.OLUnitCellR1 = R1_;
	geom.OLUnitCellR2 = R2_;
	geom.OLUnitCellAngleDeg = w_;
	geom.OLUnitCellAngleRad = w_*pi/180;
	
	% make calls to necessary update functions
	updateGeometryFields;
	if thetaPanel.showflag
		updateThetaPanel;
	end
end

function [R1, R2, w] = getGeometryUL % just for convenience
	R1 = geom.ULUnitCellR1;
	R2 = geom.ULUnitCellR2;
	w  = geom.ULUnitCellAngleDeg;
end

function [R1_, R2_, w_] = getGeometryOL
	R1_ = geom.OLUnitCellR1;
	R2_ = geom.OLUnitCellR2;
	w_  = geom.OLUnitCellAngleDeg;
end

function updateScaleRSP % this adjusts the default viewing size of the RSP
	% the resetScaleFlag indicates whether to reset the viewing range
	% the 10/-10 values here are hardcoded; change if you want a different default
	% zoom in rsp
	rsp.defaultScale = [-10 10 -10 10] * max([geom.ULUnitCellR1, geom.ULUnitCellR2, geom.OLUnitCellR1, geom.OLUnitCellR2]);
 	if rsp.resetScaleFlag
 		rsp.currentScale = rsp.defaultScale;
 	end
end

function updateScaleKSP % this adjusts the default viewing size of the KSP
	% the resetScaleFlag indicates whether to reset the viewing range
	% the 2/-2 values here are hardcoded; change if you want a different default
	% zoom in k-space
	ksp.defaultScale = [-2 2 -2 2] * max([ksp.ULk1, ksp.ULk2, ksp.OLk1, ksp.OLk2]);
 	if ksp.resetScaleFlag
 		ksp.currentScale = ksp.defaultScale;
 	end
end

function setSecondAtomPositions(UL1, UL2, OL1, OL2)
	if(any([UL1, UL2, OL1, OL2] > 1 | [UL1, UL2, OL1, OL2] < 0))
		errordlg('Use relative coordinates (0...1) for 2nd atoms!');
	else
		geom.ULSecondAtomR1 = UL1;
		geom.ULSecondAtomR2 = UL2;
		geom.OLSecondAtomR1 = OL1;
		geom.OLSecondAtomR2 = OL2;
	end
	update2ndAtomFields;
end

function setSecondAtomVisibility(boolUL, boolOL)
	rsp.ULSecondAtomIsVisible = boolUL;
	rsp.OLSecondAtomIsVisible = boolOL;
	vispn.chk1.Value = boolUL;
	vispn.chk2.Value = boolOL;
end

function setNumberOfUnitCells(nUL, nOL) % positive integers as input arguments
	rsp.ULnUnitCells = nUL;
	rsp.OLnUnitCells = nOL;
end

function setNumberOfKSpaceUnitCells(nKSpaceUL, nKSpaceOL) % positive integers as input arguments
	% 4 or 5 seems to be a good default
	% the user is never exposed to this function, but you can change this setting
	% if you want to see more/less points in ksp
	ksp.ULnUnitCells = nKSpaceUL;
	ksp.OLnUnitCells = nKSpaceOL;
end

function setRotation(TwistAngleDeg, GlobalAngleDeg)
	rot.TwistAngleDeg = TwistAngleDeg;
	rot.TwistAngleRad = TwistAngleDeg*pi/180;
	rot.GlobalAngleDeg = GlobalAngleDeg;
	rot.GlobalAngleRad = GlobalAngleDeg*pi/180;
	
	% update UI elements - unlike for many other UI elements there is no separate function for these.
	geopn.ft.String = num2str(rot.TwistAngleDeg, '%.2f');
	geopn.st.Value = rot.TwistAngleDeg;
	geopn.fg.String = num2str(rot.GlobalAngleDeg, '%.2f');
end

function recallGeometry(str) 
	rsp.resetScaleFlag = true;
	ksp.resetScaleFlag = true;
	thetaPanel.showflag = false;
	
	switch str % {'bilayer graphene (default)', 'G on h-BN', 'a-Bi on MoS2', 'a-Sb on a-Bi'}
	case 'bilayer graphene (default)'
		setGeometryUL(2.461, 2.461, 120);
		setGeometryOL(2.461, 2.461, 120);
		setRotation(10, 0);
		setSecondAtomPositions(1.0/3, 2.0/3, 1.0/3, 2.0/3);
		setSecondAtomVisibility(true, true);
		setNumberOfUnitCells (10, 10);
		setMNPQ(1, 0, 0, 1);

	% Graphene on h-BN
	case 'G on h-BN'
		setGeometryUL(2.502, 2.502, 120);
		setGeometryOL(2.461, 2.461, 120);
		setRotation(7.0, 0);
		setSecondAtomPositions(1.0/3, 2.0/3, 1.0/3, 2.0/3);
		setSecondAtomVisibility(true, true);
		setNumberOfUnitCells (12, 12);
		setMNPQ(1, 0, 0, 1);

	% a-Bi on MoS2
	case 'a-Bi on MoS2'
		setGeometryUL(3.16, 3.16, 120);
		setGeometryOL(4.53, 4.87, 90);
		setRotation(0, -90);
		setSecondAtomPositions(1.0/3, 2.0/3, 0.5, 0.44);
		setSecondAtomVisibility(false, true);
		setNumberOfUnitCells (18, 10);
		setMNPQ(1, 0, 1, 1);
		
	% a-Sb on a-Bi
	case 'a-Sb on a-Bi'
		setGeometryUL(4.54, 4.75, 90);
		setGeometryOL(4.32, 4.74, 90);
		setRotation(0, -90);
		setSecondAtomPositions(0.5, 0.47, 0.5, 0.42);
		setSecondAtomVisibility(true, true);
		setNumberOfUnitCells (10, 10);
		setMNPQ(1, 0, 1, 1);
		
	otherwise
		errordlg('Selected example not known!');
	end
	updateRspFull;
	updateKspFull;
	thetaPanel.showflag = true;
	updateThetaPanel;
end

function setMNPQ(m, n, p, q)
	thetaPanel.m = m;
	thetaPanel.n = n;
	thetaPanel.p = p;
	thetaPanel.q = q;
	updateMNPQFields;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real Space calculations - create real space lattices       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There is one function for the UL, but two for the OL. The latter ones first
% create the lattice which is unaffected by the twist (but rotated by the global
% angle), then the second function takes this existing set of points and rotates
% them additionally by the twist angle. For the UL we do it all in one function
% which is simpler to code but less fine-grained. We can afford that, however, as
% the UL is updated much less frequently (only when changing its geometry or G)!
%%%%%%%%%%%%%%%%%
function createRealSpaceUL
	G = exp(1i*rot.GlobalAngleRad);	% rotation expressed as complex phase factor
	geom.ULBaseV1 = G*geom.ULUnitCellR1;
	geom.ULBaseV2 = G*geom.ULUnitCellR2*exp(1i*geom.ULUnitCellAngleRad);	
	[ULX,ULY] = meshgrid(-rsp.ULnUnitCells : 1 : rsp.ULnUnitCells);
	rsp.ULCoords = ULX'*geom.ULBaseV1 + ULY'*geom.ULBaseV2;
end
	
function createRealSpaceGridOL % create "untwisted" OL lattice
	G = exp(1i*rot.GlobalAngleRad);	% rotation expressed as complex phase factor
	geom.untwistedOLBaseV1 = G*geom.OLUnitCellR1;
	geom.untwistedOLBaseV2 = G*geom.OLUnitCellR2*exp(1i*geom.OLUnitCellAngleRad);	
	[OLX,OLY] = meshgrid(-rsp.OLnUnitCells : 1 : rsp.OLnUnitCells);
	geom.untwistedOL = OLX'*geom.untwistedOLBaseV1 + OLY'*geom.untwistedOLBaseV2;	
	createRealSpaceLatticeOL;
end

function createRealSpaceLatticeOL % take the untwisted lattice and rotate it
	T = exp(1i*rot.TwistAngleRad);	% twist expressed as another complex phase factor
	geom.OLBaseV1 = T * geom.untwistedOLBaseV1;
	geom.OLBaseV2 = T * geom.untwistedOLBaseV2;
	rsp.OLCoords = T * geom.untwistedOL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-Space calculation - generate reciprocal bases & lattices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function generateKSpaceGridUL
	G = exp(1i*rot.GlobalAngleRad);
	
	% reciprocal "lattice constant" in Å^-1; real numbers
	ksp.ULk1 = 1/(geom.ULUnitCellR1*sin(geom.ULUnitCellAngleRad));
	ksp.ULk2 = 1/(geom.ULUnitCellR2*sin(geom.ULUnitCellAngleRad));

	% RECIPROCAL VECTORS (COMPLEX NOTATION)
	ksp.ULBaseK1 = ksp.ULk1*G*exp(1i*geom.ULUnitCellAngleRad)*(-1i); 
	ksp.ULBaseK2 = ksp.ULk2*G*1i;
	
	[ksUX,ksUY] = meshgrid(-ksp.ULnUnitCells : 1 : ksp.ULnUnitCells);
	ksp.ULCoords = ksUX'*ksp.ULBaseK1 + ksUY'*ksp.ULBaseK2;
end

function generateKSpaceGridOL
	T = exp(1i*(rot.TwistAngleRad+rot.GlobalAngleRad));
	
	% reciprocal "lattice constant" in Å^-1; real numbers
	ksp.OLk1 = 1/(geom.OLUnitCellR1*sin(geom.OLUnitCellAngleRad));
	ksp.OLk2 = 1/(geom.OLUnitCellR2*sin(geom.OLUnitCellAngleRad));

	% RECIPROCAL VECTORS (COMPLEX NOTATION)
	ksp.OLBaseK1 = ksp.OLk1*T*exp(1i*geom.OLUnitCellAngleRad)*(-1i);
	ksp.OLBaseK2 = ksp.OLk2*T*1i;
	
	[ksOX,ksOY] = meshgrid(-ksp.OLnUnitCells : 1 : ksp.OLnUnitCells);
	ksp.OLCoords = ksOX'*ksp.OLBaseK1 + ksOY'*ksp.OLBaseK2;
end

function calculateAllDiffVectors
	% ksp.U and ksp.O are trimmed from "Coords"
	ksp.U = ksp.ULCoords;
	ksp.U(abs(ksp.U) > ksp.cutoff) = iNan; % mark values larger than the cutoff
	ksp.U(ksp.ULnUnitCells+1,ksp.ULnUnitCells+1) = iNan; % mark origin (0,0)
	ksp.U = reshape(ksp.U(~isnan(ksp.U)),1,[]); % discard marked entries
	
	ksp.O = ksp.OLCoords;
	ksp.O(abs(ksp.O) > ksp.cutoff) = iNan;
	ksp.O(ksp.OLnUnitCells+1,ksp.OLnUnitCells+1) = iNan;
	ksp.O = reshape(ksp.O(~isnan(ksp.O)),1,[]);
	ksp.allDiffVectorsK = reshape(bsxfun(@minus,ksp.O.',ksp.U),1,[]);
end

function selectMPVectors
% select three smallest MP k-vectors
% conditions from the K_mnpq set given below:
%		- must not be a NaN
%		- must be Re(K)>0 (avoids redundancy of K, -K)
%		- must be smaller than k1, k2, k1_, k2_
	Kmnpq = ksp.allDiffVectorsK(real(ksp.allDiffVectorsK) >= 0);
	Kmnpq = Kmnpq(atan2d(imag(Kmnpq), real(Kmnpq)) > -90);
	Kmnpq = Kmnpq(abs(Kmnpq) <= ksp.cutoff );
	Kmnpq = sort(Kmnpq);

	if length(Kmnpq)>2
		MP.kList = [Kmnpq(1) Kmnpq(2) Kmnpq(3)];
	elseif length(Kmnpq)==2
		MP.kList = [Kmnpq(1) Kmnpq(2) iNan];
	elseif length(Kmnpq)==1
		MP.kList = [Kmnpq(1) iNan iNan];
	elseif isempty(Kmnpq) % if the length is == 0 it is easier/faster to check for an empty list!
		errordlg('No MP vectors found! Try increasing the cutoff!');
	else
		errordlg('Error in selecting the MP vectors; unhandled length of list!');
	end
	
	lambda = 1./abs(MP.kList);
	delta_deg = 180/pi*(atan(imag(MP.kList)./real(MP.kList))-rot.GlobalAngleRad)+90;
	delta_deg = mod(delta_deg+90,180)-90;
	[MP.l1, MP.l2, MP.l3] = deal(lambda(1), lambda(2), lambda(3));
	[MP.angleDeg1, MP.angleDeg2, MP.angleDeg3] = deal(delta_deg(1), delta_deg(2), delta_deg(3));
	updateMPTable;
end

function setKSPcutoff(varargin)
	if nargin == 0
		ksp.cutoff = abs(max(min([ksp.ULk1, ksp.ULk2]), min([ksp.OLk1, ksp.OLk2])))+0.001;
	elseif nargin == 1
		ksp.cutoff = abs(varargin{1});
	else
		errordlg('"setKSPcutoff" received more inputs than it can handle');
	end
	
	% update UI elements
	updateCutoffControls;
end

function setCutoffLimits % the cutoff in k, that is. The code ensures that the cutoff is always larger than the shortest reciprocal vector and can reach several higher orders
	z = max([ksp.ULnUnitCells, ksp.OLnUnitCells]); % these are not accessible from the UI but could be set manually by the user
	u = max(min(abs([ksp.ULk1, ksp.ULk2])), min(abs([ksp.OLk1, ksp.OLk2])))+0.001;
	v = z * max(abs([ksp.ULk1, ksp.ULk2, ksp.OLk1, ksp.OLk2]))+0.001;
	ksp.CutoffLimits = [u, v];
	
	% Make sure that the current position of the slider is inside the valid range!
	if MP.sk0.Value < u
		setKSPcutoff(u);
	elseif MP.sk0.Value > v
		setKSPcutoff(v);
	end
	
	% change the limits of the slider
	% NOTE: This will unfortunately affect the slider step length. :(
	MP.sk0.Min = u;
	MP.sk0.Max = v;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure and UI elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some of the UI elements feature the 'Tag' keyword for the sole purpose of
% being used as an identifier by their callback functions.
function h = addUIcomponents
	figW = 1315;
	figH = 500;
	figX = (g.ScreenSize(3)-figW)/2;
	figY = (g.ScreenSize(4)-figH)/2;
	
	% uncomment the second line, comment the first one to allow resizing of the window
	h.fig = figure('Position', [figX, figY, figW, figH] , 'Resize', 'off');
	%h.fig = figure('Position', [figX, figY, figW, figH]);

	% panel sizes:
	h.rsp = 350;		% real space
	h.ksp = 350;		% reciprocal space

	h.mainpn = uipanel('Title', 'Visualization', 'units', 'pixels',...
		'Position', [10 10 815 487], 'BorderType', 'line', 'HighlightColor', [0 0 1],...
		 'BorderWidth', 2, 'FontWeight', 'bold');
	% create panels for displaying all the point lattices later
	h.ax1 = axes(h.mainpn, 'Units', 'pixels', 'Position', [50 80 h.rsp h.rsp]);
	h.ax1i = axes(h.mainpn, 'Units', 'pixels', 'Position', [315 345 80 80]);
	h.ax2 = axes(h.mainpn, 'Units', 'pixels', 'Position', [460 80 h.ksp h.ksp]);
	
	% Togglebutton to identify (m,n), (p,q), and K_mnpq
	ksp.IdentifyKpoints = uicontrol(h.mainpn, 'Style', 'togglebutton', 'String', '<html><center>Highlight <br>M,N,P,Q',...
		'Position', [700 10 60 40], 'Callback', @identifyKpoints, 'ForegroundColor', [0 0 1], 'FontWeight', 'bold');
		
	% Reset view pushbuttons
	vispn.bResetRSP = uicontrol(h.mainpn, 'Style', 'pushbutton', 'String', 'Reset Real Space view',...
		'Position', [50 25 130 20], 'Tag', 'Reset RSP', 'Callback', @resetView);
	vispn.bResetKSP = uicontrol(h.mainpn, 'Style', 'pushbutton', 'String', 'Reset K-Space view',...
		'Position', [460 25 130 20], 'Tag', 'Reset KSP', 'Callback', @resetView);

	% geometry panel
	h.geopn = uipanel('Title','Input Geometry','units', 'pixels',...
		'Position', [835 370 470 128], 'BorderType', 'line', 'HighlightColor', [0 0.75 0],...
		 'BorderWidth', 2, 'FontWeight', 'bold');
	
	% popup menu for example selection	
	geopn.popup = uicontrol(h.geopn, 'Style', 'popup',...
		'String', {'bilayer graphene (default)', 'G on h-BN', 'a-Bi on MoS2', 'a-Sb on a-Bi'},...
		'Position', [5 90 190 20], 'Callback', @selectExample);

	% here and below: "f" in the names means "field", "s" means "slider"
	geopn.txtUL = uicontrol(h.geopn, 'Style', 'text',	'Units', 'pixels', 'Position', [5 67 95 20],...
		'String', 'Underlayer');
	geopn.txtOL = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [105 67 95 20],...
		'String', 'Overlayer');
	geopn.fr1 = uicontrol(h.geopn, 'Style', 'edit', 'Position', [45 50 50 20],...
		'String', num2str(geom.ULUnitCellR1, '%.3f'), 'Tag', 'ULR1', 'Callback', @inputGeometry);
	geopn.fr2 = uicontrol(h.geopn, 'Style', 'edit', 'Position', [45 30 50 20],...
		'String', num2str(geom.ULUnitCellR2, '%.3f'), 'Tag', 'ULR2', 'Callback', @inputGeometry);
	geopn.fw = uicontrol(h.geopn, 'Style', 'edit', 'Position', [45 10 50 20],...
		'String', num2str(geom.ULUnitCellAngleDeg, '%.2f'), 'Tag', 'ULw', 'Callback', @inputGeometry);
	geopn.txt1 = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [5 47 40 20],...
		'String', strcat('R1 (',char(197), ')'));
	geopn.txt2 = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [5 27 40 20],...
		'String', strcat('R2 (',char(197), ')'));
	geopn.txt3 = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [5 7 40 20],...
		'String', strcat(char(969), ' (', char(176), ')'));

	% overlayer
	geopn.fr1_ = uicontrol(h.geopn, 'Style', 'edit', 'Position', [145 50 50 20],...
		'String', num2str(geom.OLUnitCellR1, '%.3f'), 'Tag', 'OLR1', 'Callback', @inputGeometry);
	geopn.fr2_ = uicontrol(h.geopn, 'Style', 'edit', 'Position', [145 30 50 20],...
		'String', num2str(geom.OLUnitCellR2, '%.3f'), 'Tag', 'OLR2', 'Callback', @inputGeometry);
	geopn.fw_ = uicontrol(h.geopn, 'Style', 'edit', 'Position', [145 10 50 20],...
		'String', num2str(geom.OLUnitCellAngleDeg, '%.2f'), 'Tag', 'OLw', 'Callback', @inputGeometry);
	geopn.txt1_ = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [105 47 40 20],...
		'String', strcat('R1` (',char(197), ')'));
	geopn.txt2_ = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [105 27 40 20],...
		'String', strcat('R2` (',char(197), ')'));
	geopn.txt3_ = uicontrol(h.geopn, 'Style','text',	'Units', 'pixels', 'Position', [105 7 40 20],...
		'String', strcat(char(969), '` (', char(176), ')'));
	
	% angles - part of the geometry panel
	geopn.txt6a = uicontrol(h.geopn, 'Style', 'text', 'HorizontalAlignment', 'left',	'Units', 'pixels', 'Position', [215 25 20 20],...
		'String', strcat('0 ', char(176)));
	geopn.txt6b = uicontrol(h.geopn, 'Style', 'text', 'HorizontalAlignment', 'left',	'Units', 'pixels', 'Position', [415 25 30 20],...
		'String', strcat('360 ', char(176)));
	% twist angle theta (t)
	geopn.ft = uicontrol(h.geopn, 'Style', 'edit', 'Position', [355 75 50 20],...
		'String', num2str(rot.TwistAngleDeg, '%.2f'), 'Callback', @inputTheta);
	geopn.txt4 = uicontrol(h.geopn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels', 'Position', [265 73 90 20],...
		'String', strcat(strcat('twist angle ', char(160), char(952), '(', char(176), ')')));
	geopn.st = uicontrol(h.geopn, 'Style', 'slider', 'Min',	0, 'Max',	360, 'SliderStep', [1/3600 1/360],...
		'Value', rot.TwistAngleDeg, 'Units', 'pixels', 'Position', [205 50 240 20], 'Callback', @inputTheta);
	% global rotation angle gamma (g)
	geopn.fg = uicontrol(h.geopn, 'Style', 'edit', 'Position', [355 5 50 20],...
		'String', num2str(rot.GlobalAngleDeg, '%.2f'), 'Callback', @inputGamma);
	geopn.txt5 = uicontrol(h.geopn, 'Style','text', 'HorizontalAlignment', 'right',	'Units', 'pixels', 'Position', [215 3 140 20],...
		'String', strcat('global rotation angle',char(160), char(947), '(', char(176), ')'));
	
	% Visibility options for real space panel
	h.vispn = uipanel('Title', 'Display options', 'units', 'pixels',...
		'Position', [835 240 260 130], 'BorderType', 'line', 'HighlightColor', [0 0 1],...
		 'BorderWidth', 2, 'FontWeight', 'bold');
	
	% define number of visible unit cells (==size) for the layers
	vispn.txt1 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'left', 'Units', 'pixels',...
		'Position', [5 86 190 20], 'String', 'number of visible unit cells (*2+1)');
	vispn.txt2 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [5 69 55 20], 'String', 'underlayer');
	vispn.txt3 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [100 69 55 20], 'String', 'overlayer');
	vispn.fN = uicontrol(h.vispn, 'Style', 'edit', 'Position', [60 71 30 20],...
		'String', num2str(rsp.ULnUnitCells), 'Tag','N', 'Callback', @inputNumberOfUnitCells);
	vispn.fN_ = uicontrol(h.vispn, 'Style', 'edit', 'Position', [155 71 30 20],...
		'String', num2str(rsp.OLnUnitCells), 'Tag','N`', 'Callback', @inputNumberOfUnitCells);
	
	% checkboxes for visibility of 2nd atoms
	vispn.txt4 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'left', 'Units', 'pixels',...
		'Position', [5 45 220 20], 'String', '2nd basis atom       coordinates');
	vispn.txt5 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [5 27 50 20], 'String', 'UL show');
	vispn.txt6 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [5 5 50 20], 'String', 'OL show');
	vispn.chk1 = uicontrol(h.vispn, 'Style','checkbox', 'Units', 'pixels',...
		'Position', [55 30 20 20], 'Value', rsp.ULSecondAtomIsVisible, 'Tag', 'ULchk', 'Callback', @inputChk2ndAtom);
	vispn.chk2 = uicontrol(h.vispn, 'Style','checkbox', 'Units', 'pixels',...
		'Position', [55 10 20 20], 'Value', rsp.OLSecondAtomIsVisible, 'Tag', 'OLchk', 'Callback', @inputChk2ndAtom);
	% coordinate fields and labels
	% UL
	vispn.txt7 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [80 27 30 20], 'String', 'UL x');
	vispn.txt8 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [165 27 30 20], 'String', 'UL y');
	vispn.fXCoordUL = uicontrol(h.vispn, 'Style', 'edit', 'Position', [110 30 50 20],...
		'String', num2str(geom.ULSecondAtomR1, '%.3f'), 'Tag', 'xCooUL', 'Callback', @input2ndCoords);
	vispn.fYCoordUL = uicontrol(h.vispn, 'Style', 'edit', 'Position', [195 30 50 20],...
		'String', num2str(geom.ULSecondAtomR2, '%.3f'), 'Tag', 'yCooUL', 'Callback', @input2ndCoords);
	% OL
	vispn.txt9 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [80 5 30 20], 'String', 'OL x');
	vispn.txt10 = uicontrol(h.vispn, 'Style','text', 'HorizontalAlignment', 'center', 'Units', 'pixels',...
		'Position', [165 5 30 20], 'String', 'OL y');
	vispn.fXCoordOL = uicontrol(h.vispn, 'Style', 'edit', 'Position', [110 10 50 20],...
		'String', num2str(geom.OLSecondAtomR1, '%.3f'), 'Tag', 'xCooOL', 'Callback', @input2ndCoords);
	vispn.fYCoordOL = uicontrol(h.vispn, 'Style', 'edit', 'Position', [195 10 50 20],...
		'String', num2str(geom.OLSecondAtomR2, '%.3f'), 'Tag', 'yCooOL', 'Callback', @input2ndCoords);
	
	
	h.tablepn = uipanel('Title', 'Found moire patterns','units', 'pixels',...
		'Position', [1110 240 195 130], 'BorderType', 'line', 'HighlightColor', [1 0 0],...
		 'BorderWidth', 2, 'FontWeight', 'bold');
	% MP panel cutoff slider and field
	MP.k0txt = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [15 85 80 20], 'String', strcat('k cutoff (Å ', char(hex2dec('207B')), char(hex2dec('00B9')),')')); % this is the representation of superscript "-1"
	MP.fk0 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [95 85 50 20],...
		'String', num2str(ksp.cutoff,'%.3f'), 'Callback', @inputCutoff);
	MP.sk0 = uicontrol(h.tablepn, 'Style', 'slider', 'Min',	0.01, 'Max',	10.01, 'SliderStep', [1/100 1/10],...
		'Value',	-1, 'Units', 'pixels', 'Position', [150 5 20 105], 'Callback', @inputCutoff);

	% output: lambda, delta (3 of each)
	table.txt1 = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [5 43 17 20], 'String', 'M1');
	table.txt2 = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [5 23 17 20], 'String', 'M2');
	table.txt3 = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [5 3 17 20], 'String', 'M3');
	table.txt4 = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [30 60 35 20], 'String', strcat(char(955), ' (Å)'));
	table.txt5 = uicontrol(h.tablepn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [85 60 35 20], 'String', strcat(char(948),' (°)'));
	table.lambda1 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [25 45 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.l1,'%.2f'));
	table.lambda2 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [25 25 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.l2,'%.2f'));
	table.lambda3 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [25 5 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.l3,'%.2f'));
	table.delta1 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [85 45 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.angleDeg1,'%.2f'));
	table.delta2 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [85 25 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.angleDeg2,'%.2f'));
	table.delta3 = uicontrol(h.tablepn, 'Style', 'edit', 'Position', [85 5 50 20], 'Enable','inactive', 'BackgroundColor', [0.8 0.8 0.8],...
		'String', num2str(MP.angleDeg3,'%.2f'));

	% panel for displaying theta dependence 
	h.thetapn = uipanel('Title', 'Twist Angle Dependence','units', 'pixels',...
		'Position', [835 10 470 230], 'BorderType', 'line', 'HighlightColor', [0.9 0.8 0],...
		 'BorderWidth', 2, 'FontWeight', 'bold');

	% create axes to display data
	h.ax1pn = axes(h.thetapn,'Units', 'pixels', 'Position', [60 45 150 150]);
	h.ax2pn = axes(h.thetapn, 'Units', 'pixels', 'Position', [265 45 150 150]);

	h.pnc.fm = uicontrol(h.thetapn, 'Style', 'edit', 'Position', [435 160 30 20],...
		'String', num2str(thetaPanel.m), 'Tag', 'm', 'Callback', @inputMNPQ);
	h.pnc.fn = uicontrol(h.thetapn, 'Style', 'edit', 'Position', [435 130 30 20],...
		'String', num2str(thetaPanel.n), 'Tag', 'n', 'Callback', @inputMNPQ);
	h.pnc.fp = uicontrol(h.thetapn, 'Style', 'edit', 'Position', [435 100 30 20],...
		'String', num2str(thetaPanel.p), 'Tag', 'p', 'Callback', @inputMNPQ);
	h.pnc.fq = uicontrol(h.thetapn, 'Style', 'edit', 'Position', [435 70 30 20],...
		'String', num2str(thetaPanel.q), 'Tag', 'q', 'Callback', @inputMNPQ);
	h.pnc.txt1 = uicontrol(h.thetapn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [416 157 15 20], 'String', 'm');
	h.pnc.txt2 = uicontrol(h.thetapn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [416 127 15 20], 'String', 'n');
	h.pnc.txt3 = uicontrol(h.thetapn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [416 97 15 20], 'String', 'p');
	h.pnc.txt4 = uicontrol(h.thetapn, 'Style','text', 'HorizontalAlignment', 'right', 'Units', 'pixels',...
		'Position', [416 67 15 20], 'String', 'q');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update UI elements, handle Callback functions & inputs	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% main panels/axes - RSP and KSP
function rspUpdate
	% read and potentially adjust field of view
	if rsp.resetScaleFlag
		updateScaleRSP;
		rsp.resetScaleFlag = false;
	else
		rsp.currentScale = [h.ax1.XLim, h.ax1.YLim];
	end
	
	% check whether second atoms should be displayed
	if(rsp.ULSecondAtomIsVisible)
		secondAtomCoordinatesUL = rsp.ULCoords +...
		geom.ULSecondAtomR1*geom.ULBaseV1 +...
		geom.ULSecondAtomR2*geom.ULBaseV2;
	else
		secondAtomCoordinatesUL = [];
	end
	
	if(rsp.OLSecondAtomIsVisible)
		secondAtomCoordinatesOL = rsp.OLCoords +...
		geom.OLSecondAtomR1*geom.OLBaseV1 +...
		geom.OLSecondAtomR2*geom.OLBaseV2;
	else
		secondAtomCoordinatesOL = [];
	end
	
	% define points in order to draw the unit cells as lines
	unitCellUL = [0, geom.ULBaseV1, geom.ULBaseV1+geom.ULBaseV2, geom.ULBaseV2, 0];
	unitCellOL = [0, geom.OLBaseV1, geom.OLBaseV1+geom.OLBaseV2, geom.OLBaseV2, 0];
	
	% draw underlayer
	plot(h.ax1, real(rsp.ULCoords(:)), imag(rsp.ULCoords(:)), '.k', 'markers', rsp.markersize_1);
	hold(h.ax1, 'on');
	plot(h.ax1, real(secondAtomCoordinatesUL(:)), imag(secondAtomCoordinatesUL(:)), '.k', 'markers', rsp.markersize_1);
	plot(h.ax1, real(unitCellUL(:)), imag(unitCellUL(:)), 'k', 'LineWidth', 1.5);
	% draw overlayer
	plot(h.ax1, real(rsp.OLCoords(:)), imag(rsp.OLCoords(:)), '.r', 'markers', rsp.markersize_2);
	plot(h.ax1, real(secondAtomCoordinatesOL(:)), imag(secondAtomCoordinatesOL(:)), '.r', 'markers', rsp.markersize_2);
	plot(h.ax1, real(unitCellOL(:)), imag(unitCellOL(:)), 'r', 'LineWidth', 1.5);
	plot(h.ax1, 0, 0, 'k+', 'Markers', 0.333*ksp.markerSize1); % display a cross as the origin

	title(h.ax1, 'Real space');
	xlabel(h.ax1, 'x (Å)');
	ylabel(h.ax1, 'y (Å)');
 	axis(h.ax1, 'equal');
	
	% Implicit scale of 10x the largest lattice parameter, change if desired!
	% (The factor 10 is b/c of (1*)rsp.currentScale; see the function "updateScaleRSP")
	axis(h.ax1, rsp.currentScale);
	hold(h.ax1, 'off');
	
	% Inset panel for unit cells.
	% The numbers (-1, +4) make sure that the cell at the origin is surrounded
	% by two unit cells in each direction; if the user sets a small total number
	% of unit cells, this will fail, unfortunately...
	ul = rsp.ULnUnitCells-1 : rsp.ULnUnitCells+4;
	ol = rsp.OLnUnitCells-1 : rsp.OLnUnitCells+4;
	
	% draw the lattices around the origin and the unit cells; use smaller markers
	plot(h.ax1i, real(rsp.ULCoords(ul(:),ul(:))), imag(rsp.ULCoords(ul(:),ul(:))), '.k', 'markers', rsp.markersize_1*0.7);
	hold(h.ax1i, 'on');
	if(rsp.ULSecondAtomIsVisible)
		plot(h.ax1i, real(secondAtomCoordinatesUL(ul(1:end-1),ul(1:end-1))), imag(secondAtomCoordinatesUL(ul(1:end-1),ul(1:end-1))), '.k', 'markers', rsp.markersize_1*0.7);
	end
	plot(h.ax1i, real(unitCellUL(:)), imag(unitCellUL(:)), 'k', 'LineWidth', 1);
	plot(h.ax1i, real(rsp.OLCoords(ol(:),ol(:))), imag(rsp.OLCoords(ol(:),ol(:))), '.r', 'markers', rsp.markersize_1*0.7);
	if(rsp.OLSecondAtomIsVisible)
		plot(h.ax1i, real(secondAtomCoordinatesOL(ol(1:end-1),ol(1:end-1))), imag(secondAtomCoordinatesOL(ol(1:end-1),ol(1:end-1))), '.r', 'markers', rsp.markersize_1*0.7);
	end
	plot(h.ax1i, real(unitCellOL(:)), imag(unitCellOL(:)), 'r', 'LineWidth', 1);
	plot(h.ax1i, 0, 0, 'k+', 'Markers', 0.333*ksp.markerSize1); % display a cross as the origin

 	axis(h.ax1i, 'equal');
	% hard-coded scale of 1.5x the largest lattice parameter! change if
	% desired, see also the function "updateScaleRSP" (contributes another factor 10)
	axis(h.ax1i, 0.15*rsp.defaultScale);
	set(h.ax1i,'XTickLabel',[],'YTickLabel',[]);
	hold(h.ax1i, 'off');
end

function kspUpdate
	% read and potentially adjust field of view
	if ksp.resetScaleFlag
		updateScaleKSP;
		ksp.resetScaleFlag = false;
	else
		ksp.currentScale =	[h.ax2.XLim, h.ax2.YLim];
	end	
	
	cla(h.ax2); % this is needed as re-plotting does not delete the rectangles
	rectangle(h.ax2, 'Position', [-ksp.cutoff -ksp.cutoff 2*ksp.cutoff 2*ksp.cutoff], 'Curvature', [1 1], 'Facecolor', [0.85 0.98 0.90], 'LineStyle', 'None');
	hold(h.ax2, 'on');
	scatter(h.ax2, real(ksp.U(:)), imag(ksp.U(:)), ksp.markerSize2, 'MarkerFaceColor','k','MarkerEdgeColor','k');
	scatter(h.ax2, real(ksp.O(:)), imag(ksp.O(:)), ksp.markerSize2, 'MarkerFaceColor','r','MarkerEdgeColor','r');
	scatter(h.ax2, real(ksp.allDiffVectorsK(:)), imag(ksp.allDiffVectorsK(:)), ksp.markerSize1, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.0);
	scatter(h.ax2, real(MP.kList(:)),imag(MP.kList(:)),ksp.markerSize2, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0.0);
	plot(h.ax2, 0, 0, 'k+', 'Markers', ksp.markerSize1); % display a cross as the origin

	% show/hide the mnpq labels. Reminder: m,n are from the UL k-point; p,q from the OL k-point
	if ksp.annotationFlag
		% create the label strings
		mnstr = strcat('(', num2str(thetaPanel.m) , ',', num2str(thetaPanel.n), ')');
		pqstr = strcat('(', num2str(thetaPanel.p) , ',', num2str(thetaPanel.q), ')');
		mnpqstr = strcat('(', num2str(thetaPanel.m) , ',', num2str(thetaPanel.n), ',', num2str(thetaPanel.p), ',', num2str(thetaPanel.q), ')');
		
		% create the coordinates
		mnx = real(thetaPanel.m * ksp.ULBaseK1 + thetaPanel.n * ksp.ULBaseK2);
		mny = imag(thetaPanel.m * ksp.ULBaseK1 + thetaPanel.n * ksp.ULBaseK2);
		pqx = real(thetaPanel.p * ksp.OLBaseK1 + thetaPanel.q * ksp.OLBaseK2);
		pqy = imag(thetaPanel.p * ksp.OLBaseK1 + thetaPanel.q * ksp.OLBaseK2);
		mnpqx = pqx-mnx;
		mnpqy = pqy-mny;
		
		% draw extra markers (useful if points are not displayed due to cutoff)
		scatter(h.ax2, mnx, mny, ksp.markerSize2, 'MarkerFaceAlpha', 0.0, 'MarkerEdgeColor', 'k');
		scatter(h.ax2, pqx, pqy, ksp.markerSize2, 'MarkerFaceAlpha', 0.0, 'MarkerEdgeColor', 'r');
		scatter(h.ax2, mnpqx, mnpqy, ksp.markerSize2, 'MarkerFaceAlpha', 0.0, 'MarkerEdgeColor', 'b');
		
		% draw labels
		text(h.ax2, mnx, mny, strcat(' \leftarrow', mnstr),  'VerticalAlignment', 'baseline', 'Clipping', 'on'); % no need for 'Color', as will be black
		text(h.ax2, pqx, pqy, strcat(' \leftarrow', pqstr), 'Color', 'red', 'VerticalAlignment', 'baseline', 'Clipping', 'on');
		text(h.ax2, mnpqx, mnpqy, strcat(' \leftarrow', mnpqstr), 'Color', 'blue', 'VerticalAlignment', 'baseline', 'Clipping', 'on');	
	end
	
	title(h.ax2, 'Reciprocal space');
	xlabel(h.ax2, 'k_x (Å^-^1)')
	ylabel(h.ax2, 'k_y (Å^-^1)')
	box(h.ax2, 'on');
 	axis(h.ax2, 'equal');
	axis(h.ax2, ksp.currentScale);
	hold(h.ax2, 'off');
	drawnow;
end

% theta-dependence panel
function updateThetaPanel
	t = -180:0.1:180; % theta range, step size 0.1 degree
	%t = t(t~=0); % removing the zero value as it can create numerical discontinuities and singularities
	R1  = geom.ULUnitCellR1;
	R2  = geom.ULUnitCellR2;
	R1_ = geom.OLUnitCellR1;
	R2_ = geom.OLUnitCellR2;
	omega  = geom.ULUnitCellAngleDeg;
	omega_ = geom.OLUnitCellAngleDeg;
	mnpq = [thetaPanel.m, thetaPanel.n, thetaPanel.p, thetaPanel.q];

	% discardZeroFlag (last arg in the function call), by default set to "false"
	[ lambda, delta ] = twistAngleDependence( R1, R2, omega, R1_, R2_ , omega_ , mnpq, t, false );
	
	% display
	plot(h.ax1pn, t, lambda,'.', 'MarkerSize', 3);
	plot(h.ax2pn, t, delta,'.', 'MarkerSize', 3);
	xlabel(h.ax1pn, '\theta (°)');
	ylabel(h.ax1pn, '\lambda(Å)');
	title(h.ax1pn, '\lambda(\theta)');
	xlim(h.ax1pn, [-185, 185]);
	h.ax1pn.XTick = -180:60:180;

	xlabel(h.ax2pn, '\theta (°)');
	ylabel(h.ax2pn, '\delta(°)');
	title(h.ax2pn, '\delta(\theta)');
	axis(h.ax2pn, [-185, 185, -91, 91]);
	h.ax2pn.XTick = -180:60:180;
	h.ax2pn.YTick = -90:30:90;
end

function updateMNPQFields
	h.pnc.fm.String = num2str(thetaPanel.m);
	h.pnc.fn.String = num2str(thetaPanel.n);
	h.pnc.fp.String = num2str(thetaPanel.p);
	h.pnc.fq.String = num2str(thetaPanel.q);
end

% update input fields/sliders of main figure
function updateMPTable
	table.lambda1.String = num2str(MP.l1,'%.2f');
	table.lambda2.String = num2str(MP.l2,'%.2f');
	table.lambda3.String = num2str(MP.l3,'%.2f');
	table.delta1.String = num2str(MP.angleDeg1,'%.2f');
	table.delta2.String = num2str(MP.angleDeg2,'%.2f');
	table.delta3.String = num2str(MP.angleDeg3,'%.2f');
end

function updateGeometryFields
	geopn.fr1.String = num2str(geom.ULUnitCellR1, '%.3f');
	geopn.fr2.String = num2str(geom.ULUnitCellR2, '%.3f');
	geopn.fw.String = num2str(geom.ULUnitCellAngleDeg, '%.2f');
	geopn.fr1_.String = num2str(geom.OLUnitCellR1, '%.3f');
	geopn.fr2_.String = num2str(geom.OLUnitCellR2, '%.3f');
	geopn.fw_.String = num2str(geom.OLUnitCellAngleDeg, '%.2f');
end

function update2ndAtomFields
	vispn.fXCoordUL.String = num2str(geom.ULSecondAtomR1, '%.3f');
	vispn.fYCoordUL.String = num2str(geom.ULSecondAtomR2, '%.3f');
	vispn.fXCoordOL.String = num2str(geom.OLSecondAtomR1, '%.3f');
	vispn.fYCoordOL.String = num2str(geom.OLSecondAtomR2, '%.3f');
end

function updateNumOfUnitCellFields
	vispn.fN.String = num2str(rsp.ULnUnitCells);
	vispn.fN_.String = num2str(rsp.OLnUnitCells);
end

function updateCutoffControls
	MP.fk0.String = num2str(ksp.cutoff, '%.2f');
	MP.sk0.Value = ksp.cutoff;
end

% there is no function for theta and gamma since they both rely on the same
% function to set the rotations and it would be pointless to create a separate
% function for only gamma to set a single UI element. Instead this is done in the
% setRotation function.

% fine-grained different variants of what to do when updating the main panels
function updateRspFull % triggers a FULL real-space update
	createRealSpaceUL;
	createRealSpaceGridOL;
	createRealSpaceLatticeOL;
	rspUpdate;
	%drawnow;
end

function updateRspUL % triggers a partial real-space update for the UL
	createRealSpaceUL;
	rspUpdate;
	%drawnow;
end

function updateRspOL % triggers a partial real-space update for the OL without twist
	createRealSpaceGridOL;
	createRealSpaceLatticeOL;
	rspUpdate;
	%drawnow;
end

function updateRspOLTwist % triggers a partial real-space update for the OL doing ONLY the twist
	createRealSpaceLatticeOL;
	rspUpdate;
	%drawnow;
end

function updateKspFull(varargin) % Triggers a FULL k-space update.
	% Providing any arbitrary input argument causes the cutoff limit setting to be
	% skipped.
	generateKSpaceGridUL;
	generateKSpaceGridOL;
	if nargin == 0
		setCutoffLimits;
	end
	calculateAllDiffVectors;
	selectMPVectors;
	kspUpdate;
	%drawnow;
end

function updateKspUL % triggers a partial k-space update for the UL
	generateKSpaceGridUL;
	setCutoffLimits;
	calculateAllDiffVectors;
	selectMPVectors;
	kspUpdate;
	%drawnow;
end

function updateKspOL % triggers a partial k-space update for the OL
	generateKSpaceGridOL;
	setCutoffLimits;
	calculateAllDiffVectors;
	selectMPVectors;
	kspUpdate;
	%drawnow;
end

% input handling
function selectExample(source, ~) % handle input/callback from dropdown menu
	recallGeometry(source.String{source.Value});
	updateNumOfUnitCellFields;
end

function inputCutoff(source, ~) % handles input from k_cutoff field & slider
	if strcmp(source.Style, 'slider')
		newCutoffValue = source.Value;
	elseif strcmp(source.Style, 'edit')
		newCutoffValue = str2double(source.String);
	else
		errordlg('Function "inputCutoff" called from wrong source!')
	end

	limitsMessage = 'Limit reached! If you want other cutoff values, please edit the function "setCutoffLimits"';
	if newCutoffValue < MP.sk0.Min
		newCutoffValue = MP.sk0.Min;
		warndlg(limitsMessage);
	elseif newCutoffValue > MP.sk0.Max
		newCutoffValue = MP.sk0.Max;
		warndlg(limitsMessage);
	end	
	
	setKSPcutoff(newCutoffValue);
	updateKspFull('skip setting cutoff limits'); % see definition of updateKspFull
end

function inputTheta(source, ~) % handles input from theta field & slider
	if strcmp(source.Style, 'slider')
		newThetaValueDeg = source.Value;
	elseif strcmp(source.Style, 'edit')
		newThetaValueDeg = str2double(source.String);
	else
		warndlg('Wrong source when calling "inputTheta"!');
	end
	
	% limit theta for the slider
	if newThetaValueDeg < geopn.st.Min
		newThetaValueDeg = geopn.st.Min;
	elseif newThetaValueDeg > geopn.st.Max
		newThetaValueDeg = geopn.st.Max;
	end
	
	% set the new twist angle value
	setRotation(newThetaValueDeg, rot.GlobalAngleDeg);

	%update the data structures and displays
	updateRspOLTwist;
	updateKspOL;
end

function inputGamma(source, ~) % handles input from the gamma field
	setRotation(rot.TwistAngleDeg, str2double(source.String));
	% now update the full data structure (i.e. both layers) and displays
	updateRspFull;
	updateKspFull;
end

function inputGeometry(source, ~) % handle input for the unit cell geometries
switch source.Tag
	case 'ULR1'
		[~, R2, w] = getGeometryUL;
		R1 = str2double(source.String);
		setGeometryUL(R1, R2, w);
		updateRspUL;
		updateKspUL;
	case 'ULR2'
		[R1, ~, w] = getGeometryUL;
		R2 = str2double(source.String);
		setGeometryUL(R1, R2, w);
		updateRspUL;
		updateKspUL;
	case 'ULw'
		[R1, R2, ~] = getGeometryUL;
		w = str2double(source.String);
		setGeometryUL(R1, R2, w);
		updateRspUL;
		updateKspUL;
	case 'OLR1'
		[~, R2, w] = getGeometryOL;
		R1 = str2double(source.String);
		setGeometryOL(R1, R2, w);
		updateRspOL;
		updateKspOL;
	case 'OLR2'
		[R1, ~, w] = getGeometryOL;
		R2 = str2double(source.String);
		setGeometryOL(R1, R2, w);
		updateRspOL;
		updateKspOL;
	case 'OLw'
		[R1, R2, ~] = getGeometryOL;
		w = str2double(source.String);
		setGeometryOL(R1, R2, w);
		updateRspOL;
		updateKspOL;
	otherwise
		errordlg('Wrong source to "inputGeometry": %s', source.Tag); % might be useful when adding more examples
	end
end

function inputNumberOfUnitCells(source, ~) % handle input from the number of unit cell fiels
	n = str2double(source.String);
	if n<4
		errordlg('Please enter a number >= 4!');
	elseif n>100
		errordlg('Warning! You are attempting to display a huge number of unit cells. Please enter a smaller number <= 100!');
	else

		nUL = rsp.ULnUnitCells;
		nOL = rsp.OLnUnitCells;
		switch source.Tag
			case 'N' % note the subtle difference between N and N` (the "`")
				nUL = n;
			case 'N`'
				nOL = n;
			otherwise
		end
		setNumberOfUnitCells(nUL, nOL);
		updateNumOfUnitCellFields;
		updateRspFull;
	end
end

function inputChk2ndAtom(source, ~) % checkboxes for visibility of 2nd atoms
	boolUL = rsp.ULSecondAtomIsVisible;
	boolOL = rsp.OLSecondAtomIsVisible;
	switch source.Tag
		case 'ULchk'
			boolUL = source.Value;
		case 'OLchk'
			boolOL = source.Value;
		otherwise
			errordlg('Bad usage of 2nd atom checkbox function!');
	end
	setSecondAtomVisibility(boolUL, boolOL);
	rspUpdate;
end

function input2ndCoords(source, ~) % handle input for 2nd atom coordinates
	UL1 = geom.ULSecondAtomR1;
	UL2 = geom.ULSecondAtomR2;
	OL1 = geom.OLSecondAtomR1;
	OL2 = geom.OLSecondAtomR2;
	
	switch source.Tag
		case 'xCooUL'
			UL1= str2double(source.String);
		case 'yCooUL'
			UL2= str2double(source.String);
		case 'xCooOL'
			OL1= str2double(source.String);
		case 'yCooOL'
			OL2= str2double(source.String);
		otherwise
			errordlg('Bad usage of 2nd coordinate input function!');
	end
	setSecondAtomPositions(UL1, UL2, OL1, OL2);
	rspUpdate;
end 

function identifyKpoints(source, event) % handle button press to show/hide Kmnpq in KSP
	ksp.annotationFlag = source.Value;
	kspUpdate;
end

function inputMNPQ(source, ~) % handle input from the mnpq fields
	m = thetaPanel.m;
	n = thetaPanel.n;
	p = thetaPanel.p;
	q = thetaPanel.q;
	
	switch source.Tag
		case 'm'
			m = str2double(source.String);
		case 'n'
			n = str2double(source.String);
		case 'p'
			p = str2double(source.String);
		case 'q'
			q = str2double(source.String);
		otherwise
			errordlg('MNPQ input from wrong source!');
	end
	setMNPQ(m, n, p, q);
	updateMNPQFields;
	updateThetaPanel;
	if ksp.annotationFlag
		kspUpdate;
	end
end 

function resetView(source, ~) % handle the "reset view" buttons for RSP & KSP
	switch source.Tag
		case 'Reset RSP'
			rsp.resetScaleFlag = true;
			rspUpdate;
		case 'Reset KSP'
			ksp.resetScaleFlag = true;
			kspUpdate;
	end
end

function myweb(~,~)
	myurl = 'https://doi.org/10.1088/2053-1583/ab5470';
	web(myurl);
end
end