function [T, Y] = Rosenbrock_nonAdaptive(Function, Tspan, Y0, Options, ICNTRL)

% Parse ODE options
Jacobian = odeget(Options, 'Jacobian');
if isempty(Jacobian)
    error('A Jacobian function is required.');
end

% Get problem size
steps = length(Tspan);
N = length(Y0);

% Integrate
Y = zeros(steps,N);
T = zeros(steps,1);
Y(1,:) = Y0;
T(1) = Tspan(1);

for i=2:steps
    [T(i), Y(i,:)] = RosenbrockIntegrator(N, Y(i-1,:), T(i-1),Tspan(i), ...
            Function, Jacobian, ICNTRL);
end
return

function [ T, Y,IERR,ISTATUS,RSTATUS] = RosenbrockIntegrator(N, Y, Tstart, Tend, ...
    OdeFunction, OdeJacobian, ICNTRL)
% Advances the ODE system defined by OdeFunction, OdeJacobian and Y
% from time=Tstart to time=Tend


%~~~>   Parameters
global ZERO ONE DeltaMin
ZERO = 0.0; ONE = 1.0; DeltaMin = 1.0E-14; % DeltaMin = 1.0E-5

%~~~>  Autonomous or time dependent ODE. Default is time dependent.
Autonomous = ~(ICNTRL(1) == 0);

%~~~>   Initialize the particular Rosenbrock method selected
switch (ICNTRL(3))
    case (1)
        ros_Param = Ros2;
    case (2)
        ros_Param = Ros3;
    case (3)
        ros_Param = Ros4;
    case {0,4}
        ros_Param = Rodas3;
    case (5)
        ros_Param = Rodas4;
    otherwise
        disp(['Unknown Rosenbrock method: ICNTRL(3)=',num2str(ICNTRL(3))]);
        IERR = ros_ErrorMsg(-2,Tstart,ZERO);
        return
end
ros_S     = ros_Param{1};
rosMethod = ros_Param{2};
ros_A     = ros_Param{3};
ros_C     = ros_Param{4};
ros_M     = ros_Param{5};
ros_E     = ros_Param{6};
ros_Alpha = ros_Param{7};
ros_Gamma = ros_Param{8};
ros_NewF  = ros_Param{10};
ros_Name  = ros_Param{11};

%~~~>  Unit roundoff (1+Roundoff>1)
Roundoff = eps;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Template for the implementation of a generic Rosenbrock method
%      defined by ros_S (no of stages)
%      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~ Local variables
Ynew  = zeros(N,1);
Fcn0  = zeros(N,1);
Fcn   = zeros(N,1);
K     = zeros(N,ros_S);
dFdT  = zeros(N,1);
Jac0  = sparse(N,N);
Ghimj = sparse(N,N);
Yerr  = zeros(N,1);
Pivot = zeros(N,1);

% Write H here
%~~~>  Initial preparations
T = Tstart;
H = Tend - Tstart;

if Tend >= Tstart
    Direction = +1;
else
    Direction = -1;
end
H = Direction*H;

%~~~>   Compute the function at current time
Fcn0 = OdeFunction(T, Y);

%~~~>  Compute the function derivative with respect to T
if ~Autonomous
	dFdT = ros_FunTimeDerivative (T, OdeFunction, Roundoff, Y, Fcn0);
end

%~~~>   Compute the Jacobian at current time
Jac0 = OdeJacobian(T,Y);

%~~~>  Repeat step calculation until current step accepted
    
[H, Ghimj, Pivot, Singular] = ros_PrepareMatrix(N, H, Direction, ros_Gamma(1), Jac0);

% For the 1st istage the function has been computed previously
Fcn = Fcn0;
K(:,1) = Fcn;
if (~Autonomous) && (ros_Gamma(1) ~= ZERO)
	HG = Direction*H*ros_Gamma(1);
	K(:,1) = K(:,1) + HG * dFdT;
end
K(:,1) = ros_Solve(Ghimj, Pivot, K(:,1));
        
%~~~>   Compute the remaining stages
for istage=2:ros_S

% istage>1 and a new function evaluation is needed
%if ros_NewF(istage)
%	Ynew = Y;
%             for j=1:istage-1
%                  Ynew = Ynew + ros_A((istage-1)*(istage-2)/2+j)*K(:,j)';
%             end
%        Tau = T + ros_Alpha(istage)*Direction*H;
%        Fcn = OdeFunction(Tau,Ynew);
%end

K(:,istage) = Fcn;
	for j=1:istage-1
        	HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H);
                K(:,istage) = K(:,istage) + HC * K(:,j);
        end
        if (~Autonomous) && (ros_Gamma(istage) ~= ZERO)
                HG = Direction*H*ros_Gamma(istage);
                K(:,istage) = K(:,istage) + HG * dFdT;
        end

        % Linear system is solved here with MATLAB '\' operator
        % instead of LU decomposition.
        K(:,istage) = ros_Solve(Ghimj, Pivot, K(:,istage));
end

        %~~~>  Compute the new solution
        Ynew = Y;
        for j=1:ros_S
            Ynew = Ynew + ros_M(j) * K(:,j)';
        end
        T = T + Direction*H;
        Y = Ynew;
        Yerr = zeros(N,1);
        IERR = 1;  %~~~> The integration was successful
return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ dFdT ] = ros_FunTimeDerivative (T, OdeFunction, Roundoff, Y, Fcn0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~> The time partial derivative of the function by finite differences
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~>  Statistics on the work performed by the Rosenbrock method
global Nfun ISTATUS

%~~~>   Parameters
global ZERO ONE DeltaMin

Delta = sqrt(Roundoff) * max(DeltaMin, abs(T));
dFdT = OdeFunction(T+Delta,Y);
dFdT = (dFdT - Fcn0) / Delta;
return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [H, Ghimj, Pivot, Singular] = ...
    ros_PrepareMatrix(N, H, Direction, gam, Jac0)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Prepares the LHS matrix for stage calculations
%  1.  Construct Ghimj = 1/(H*ham) - Jac0
%      "(Gamma H) Inverse Minus Jacobian"
%  2.  LU decomposition not performed here because
%      MATLAB solves the linear system with '\'.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ghimj = -Jac0;
ghinv = 1.0/(Direction*H*gam);
for i=1:N
    Ghimj(i,i) = Ghimj(i,i)+ghinv;
end

Pivot(1) = 1;
Singular = false;

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [X] = ros_Solve(JVS, Pivot, X)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Template for the forward/backward substitution 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global Nsol

X = JVS\X;
Pivot(1) = 1;
return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros2()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- AN L-STABLE METHOD, 2 stages, order 2
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

g = 1.0 + 1.0/sqrt(2.0);
rosMethod = 1;
%~~~> Name of the method
ros_Name = 'ROS-2';
%~~~> Number of stages
ros_S = 2;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = (1.0)/g;
ros_C(1) = (-2.0)/g;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
%~~~> M_i = Coefficients for new step solution
ros_M(1)= (3.0)/(2.0*g);
ros_M(2)= (1.0)/(2.0*g);
% E_i = Coefficients for error estimator
ros_E(1) = 1.0/(2.0*g);
ros_E(2) = 1.0/(2.0*g);
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus one
ros_ELO = 2.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 1.0;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = g;
ros_Gamma(2) =-g;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros3()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 2;
%~~~> Name of the method
ros_Name = 'ROS-3';
%~~~> Number of stages
ros_S = 3;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1)= 1.0;
ros_A(2)= 1.0;
ros_A(3)= 0.0;

ros_C(1) = -0.10156171083877702091975600115545E+01;
ros_C(2) =  0.40759956452537699824805835358067E+01;
ros_C(3) =  0.92076794298330791242156818474003E+01;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
ros_NewF(3) = false;
%~~~> M_i = Coefficients for new step solution
ros_M(1) =  0.1E+01;
ros_M(2) =  0.61697947043828245592553615689730E+01;
ros_M(3) = -0.42772256543218573326238373806514;
% E_i = Coefficients for error estimator
ros_E(1) =  0.5;
ros_E(2) = -0.29079558716805469821718236208017E+01;
ros_E(3) =  0.22354069897811569627360909276199;
%~~~> ros_ELO = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO = 3.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1)= 0.0;
ros_Alpha(2)= 0.43586652150845899941601945119356;
ros_Alpha(3)= 0.43586652150845899941601945119356;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1)= 0.43586652150845899941601945119356;
ros_Gamma(2)= 0.24291996454816804366592249683314;
ros_Gamma(3)= 0.21851380027664058511513169485832E+01;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Ros4()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
%     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
%
%      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
%      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
%      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
%      SPRINGER-VERLAG (1990)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 3;
%~~~> Name of the method
ros_Name = 'ROS-4';
%~~~> Number of stages
ros_S = 4;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.2000000000000000E+01;
ros_A(2) = 0.1867943637803922E+01;
ros_A(3) = 0.2344449711399156;
ros_A(4) = ros_A(2);
ros_A(5) = ros_A(3);
ros_A(6) = 0.0;

ros_C(1) =-0.7137615036412310E+01;
ros_C(2) = 0.2580708087951457E+01;
ros_C(3) = 0.6515950076447975;
ros_C(4) =-0.2137148994382534E+01;
ros_C(5) =-0.3214669691237626;
ros_C(6) =-0.6949742501781779;
%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1)  = true;
ros_NewF(2)  = true;
ros_NewF(3)  = true;
ros_NewF(4)  = false;
%~~~> M_i = Coefficients for new step solution
ros_M(1) = 0.2255570073418735E+01;
ros_M(2) = 0.2870493262186792;
ros_M(3) = 0.4353179431840180;
ros_M(4) = 0.1093502252409163E+01;
%~~~> E_i  = Coefficients for error estimator
ros_E(1) =-0.2815431932141155;
ros_E(2) =-0.7276199124938920E-01;
ros_E(3) =-0.1082196201495311;
ros_E(4) =-0.1093502252409163E+01;
%~~~> ros_ELO  = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO  = 4.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 0.1145640000000000E+01;
ros_Alpha(3) = 0.6552168638155900;
ros_Alpha(4) = ros_Alpha(3);
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.5728200000000000;
ros_Gamma(2) =-0.1769193891319233E+01;
ros_Gamma(3) = 0.7592633437920482;
ros_Gamma(4) =-0.1049021087100450;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Rodas3()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% --- A STIFFLY-STABLE METHOD, 4 stages, order 3
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 4;
%~~~> Name of the method
ros_Name = 'RODAS-3';
%~~~> Number of stages
ros_S = 4;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:
%       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.0;
ros_A(2) = 2.0;
ros_A(3) = 0.0;
ros_A(4) = 2.0;
ros_A(5) = 0.0;
ros_A(6) = 1.0;

ros_C(1) = 4.0;
ros_C(2) = 1.0;
ros_C(3) =-1.0;
ros_C(4) = 1.0;
ros_C(5) =-1.0;
ros_C(6) =-(8.0/3.0);

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1)  = true;
ros_NewF(2)  = false;
ros_NewF(3)  = true;
ros_NewF(4)  = true;
%~~~> M_i = Coefficients for new step solution
ros_M(1) = 2.0;
ros_M(2) = 0.0;
ros_M(3) = 1.0;
ros_M(4) = 1.0;
%~~~> E_i  = Coefficients for error estimator
ros_E(1) = 0.0;
ros_E(2) = 0.0;
ros_E(3) = 0.0;
ros_E(4) = 1.0;
%~~~> ros_ELO  = estimator of local order - the minimum between the
%    main and the embedded scheme orders plus 1
ros_ELO  = 3.0;
%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.0;
ros_Alpha(2) = 0.0;
ros_Alpha(3) = 1.0;
ros_Alpha(4) = 1.0;
%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.5;
ros_Gamma(2) = 1.5;
ros_Gamma(3) = 0.0;
ros_Gamma(4) = 0.0;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [ params ] = Rodas4()
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
%
%      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
%      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
%      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
%      SPRINGER-VERLAG (1996)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rosMethod = 5;
%~~~> Name of the method
ros_Name = 'RODAS-4';
%~~~> Number of stages
ros_S = 6;

%~~~> Y_stage_i ~ Y( T + H*Alpha_i )
ros_Alpha(1) = 0.000;
ros_Alpha(2) = 0.386;
ros_Alpha(3) = 0.210;
ros_Alpha(4) = 0.630;
ros_Alpha(5) = 1.000;
ros_Alpha(6) = 1.000;

%~~~> Gamma_i = \sum_j  gamma_{i,j}
ros_Gamma(1) = 0.2500000000000000;
ros_Gamma(2) =-0.1043000000000000;
ros_Gamma(3) = 0.1035000000000000;
ros_Gamma(4) =-0.3620000000000023E-01;
ros_Gamma(5) = 0.0;
ros_Gamma(6) = 0.0;

%~~~> The coefficient matrices A and C are strictly lower triangular.
%   The lower triangular (subdiagonal) elements are stored in row-wise order:
%   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
%   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
%                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

ros_A(1) = 0.1544000000000000E+01;
ros_A(2) = 0.9466785280815826;
ros_A(3) = 0.2557011698983284;
ros_A(4) = 0.3314825187068521E+01;
ros_A(5) = 0.2896124015972201E+01;
ros_A(6) = 0.9986419139977817;
ros_A(7) = 0.1221224509226641E+01;
ros_A(8) = 0.6019134481288629E+01;
ros_A(9) = 0.1253708332932087E+02;
ros_A(10) =-0.6878860361058950;
ros_A(11) = ros_A(7);
ros_A(12) = ros_A(8);
ros_A(13) = ros_A(9);
ros_A(14) = ros_A(10);
ros_A(15) = 1.0;

ros_C(1) =-0.5668800000000000E+01;
ros_C(2) =-0.2430093356833875E+01;
ros_C(3) =-0.2063599157091915;
ros_C(4) =-0.1073529058151375;
ros_C(5) =-0.9594562251023355E+01;
ros_C(6) =-0.2047028614809616E+02;
ros_C(7) = 0.7496443313967647E+01;
ros_C(8) =-0.1024680431464352E+02;
ros_C(9) =-0.3399990352819905E+02;
ros_C(10) = 0.1170890893206160E+02;
ros_C(11) = 0.8083246795921522E+01;
ros_C(12) =-0.7981132988064893E+01;
ros_C(13) =-0.3152159432874371E+02;
ros_C(14) = 0.1631930543123136E+02;
ros_C(15) =-0.6058818238834054E+01;

%~~~> M_i = Coefficients for new step solution
ros_M(1) = ros_A(7);
ros_M(2) = ros_A(8);
ros_M(3) = ros_A(9);
ros_M(4) = ros_A(10);
ros_M(5) = 1.0;
ros_M(6) = 1.0;

%~~~> E_i  = Coefficients for error estimator
ros_E(1) = 0.0;
ros_E(2) = 0.0;
ros_E(3) = 0.0;
ros_E(4) = 0.0;
ros_E(5) = 0.0;
ros_E(6) = 1.0;

%~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
%   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
ros_NewF(1) = true;
ros_NewF(2) = true;
ros_NewF(3) = true;
ros_NewF(4) = true;
ros_NewF(5) = true;
ros_NewF(6) = true;

%~~~> ros_ELO  = estimator of local order - the minimum between the
%        main and the embedded scheme orders plus 1
ros_ELO = 4.0;

params = { ros_S, rosMethod, ros_A, ros_C, ros_M, ros_E, ...
    ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, ros_Name };

return

% End of INTEGRATE function
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


