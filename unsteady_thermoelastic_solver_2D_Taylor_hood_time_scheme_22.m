function [u1h,u2h,temp_h] = unsteady_thermoelastic_solver_2D_Taylor_hood_time_scheme_22(basis_type_temp, basis_type_disp, Gauss_point_number_1D,Gauss_point_number_2D,end_time,start_time,dt,theta,shape,refine_times)
%Wedam Nyaaba 10/30/2016

% result = unsteady_thermoelastic_solver_2D_Taylor_hood_time_scheme_22(left, right,bottom,top,h_partition, basis_type, Gauss_point_number, ... ) solves the 2D unsteady thermoelastic equation presented in the notes with an assumed heat source (H)
% Two (FE)basis types: Linear and Quadratic cases. basis_type = 1 is the Linear case; basis_type = 2 represents the quadratic case.
% left: the left-most boundary node of the domain. 
% right: right-most boundary node of the domain.
% h is the partition size of the mesh.
% Gauss_point_number: the type of Gauss nodes used in the Gauss quadrature approximation of the basis functions and coefficient functions.
% Pb matrix: stores the coordinates of the jth FE nodes.
% Tb matrix: stores the global node indices of the mesh nodes of the nth FE mesh element.
% P matrix: stores the jth coordinates of the nth mesh elements
% T matrix: stores the global node indices of the mesh nodes of the nth mesh element

% N1_basis,N2_basis: the N1 and N2 for the FE basis functions,not the partition.
% N1: the number of the sub-intervals of the partition in x-direction.
% N2: the number of the sub-intervals of the partition in y-direction.
% N_elements: the number of the mesh partitions. For 1D Linear FE, N_basis = N_elements
% number_of_lb: the number of local basis functions at the element level.
% function_lambda:function of Lames first constant
% function_mu: function of Lame second constant
% function_alpha: function of Lame constants coefficient of thermal expansion
% function_lambda_alpha: function of the product of lambda and alpha
% function_mu_alpha: function of the product of mu and alpha
% function_K: function of isotropic heat conductivity
% function_density_spec_ht: function of the product of mass density and specific heat capacity 
% function_density : function of mass density
% function_initial_conditions_disp_1: function 

% function_f is the coefficient function (f) portion of the vector integrand
%function_g_Dirichlet: the Dirichelet boundary function in my notes "Notes for tool box of standard triangular FE" section 1-1.
%function_p_tilde is: c*p in the given notes
% h_partition is a vector of [h1 h2]
% N_t is the number of time steps
% fixed_matrix_1 and fixed_matrix_2 are the fixed (chosen because coefficients are time-independent) matrix coefficients of X^m+1 and X^m, respectively

% h1 = h_partition(1);
% h2 = h_partition(2);
% N1 = (right-left)/h1;
% N2 = (top-bottom)/h2;

[P,T,Pb,Tb,boundary_nodes_mesh,boundary_nodes_basis,boundary_edges_heat,boundary_edges_stress1]=GenerateMeshData_withRobinBc(shape,1,refine_times);

N_elements = size(Tb,2);
%% Define FE basis type parameters
if basis_type_temp==2 
%     N1_basis_temp = 2*N1;
%     N2_basis_temp = 2*N2;
%   N_basis_temp = (N1_basis_temp+1)*(N2_basis_temp+1);
  number_of_lb_temp = 6;
  basis_type_trial_temp = 2;
  basis_type_test_temp = 2;
 [P,T,Pb_temp,Tb_trial_temp,boundary_nodes_mesh_temp,boundary_nodes_basis_heat,boundary_edges_heat,boundary_edges_stress]=GenerateMeshData_withRobinBc(shape,basis_type_temp,refine_times);
  N_basis_temp = size(Pb_temp,2);
  Tb_test_temp = Tb_trial_temp;
   
elseif basis_type_temp==1 
%     N1_basis_temp = N1;
%     N2_basis_temp = N2;
%   N_basis_temp = (N1_basis_temp+1)*(N2_basis_temp+1);
  number_of_lb_temp = 3;
  basis_type_trial_temp = 1;
  basis_type_test_temp = 1;
 [P,T,Pb_temp,Tb_trial_temp,boundary_nodes_mesh_temp,boundary_nodes_basis_heat,boundary_edges_heat,boundary_edges_stress]=GenerateMeshData_withRobinBc(shape,basis_type_temp,refine_times);
 N_basis_temp = size(Pb_temp,2);
 Tb_test_temp = Tb_trial_temp;
    
end
% Defining quadratic FE parameters for displacement component
if  basis_type_disp==2
%   N1_basis_disp = 2*N1;
%   N2_basis_disp = 2*N2;
%   N_basis_disp = (N1_basis_disp+1)*(N2_basis_disp+1);             
  number_of_lb_disp = 6;
  basis_type_trial_disp= 2;
  basis_type_test_disp= 2;
  [P,T,Pb_disp,Tb_trial_disp,boundary_nodes_mesh_disp,boundary_nodes_basis_disp,boundary_edges_heat,boundary_edges_stress]=GenerateMeshData_withRobinBc(shape,basis_type_disp,refine_times);
    N_basis_disp = size(Pb_disp,2);
  Tb_test_disp = Tb_trial_disp;
elseif  basis_type_disp==1
%   N1_basis_disp = N1;
%   N2_basis_disp = N2;
%   N_basis_disp = (N1_basis_disp+1)*(N2_basis_disp+1);             
  number_of_lb_disp = 3;
  basis_type_trial_disp= 1;
  basis_type_test_disp= 1;
 [P,T,Pb_disp,Tb_trial_disp,boundary_nodes_mesh_disp,boundary_nodes_basis_disp,boundary_edges_heat,boundary_edges_stress]=GenerateMeshData_withRobinBc(shape,basis_type_disp,refine_times);
 N_basis_disp = size(Pb_disp,2);
 Tb_test_disp = Tb_trial_disp;
end
% Obtain Gauss nodes

[ref_Gauss_weight_1D,ref_Gauss_nodes_1D]=generate_Gauss_reference_1D(Gauss_point_number_1D);
[ref_Gauss_weight,ref_Gauss_nodes]=generate_Gauss_reference_2D_triangle(Gauss_point_number_2D);

%% Assembly Stiffness Matrix

A1 = assemble_matrix_volume_int_2D('function_lambda',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,1,0,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A2 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,1,0,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A3 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,1,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A4 = assemble_matrix_volume_int_2D('function_lambda',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,1,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A5 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,1,0,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A6 = assemble_matrix_volume_int_2D('function_lambda_alpha',P,T,Tb_trial_temp,Tb_test_disp,number_of_lb_temp,number_of_lb_disp,N_basis_temp,N_basis_disp,N_elements,basis_type_trial_temp,0,0,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A7 = assemble_matrix_volume_int_2D('function_mu_alpha',P,T,Tb_trial_temp,Tb_test_disp,number_of_lb_temp,number_of_lb_disp,N_basis_temp,N_basis_disp,N_elements,basis_type_trial_temp,0,0,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A8 = assemble_matrix_volume_int_2D('function_lambda',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,1,0,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A9 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,1,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A10 = assemble_matrix_volume_int_2D('function_lambda',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,1,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A11 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,1,0,basis_type_test_disp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A12 = assemble_matrix_volume_int_2D('function_mu',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,1,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A13 = assemble_matrix_volume_int_2D('function_lambda_alpha',P,T,Tb_trial_temp,Tb_test_disp,number_of_lb_temp,number_of_lb_disp,N_basis_temp,N_basis_disp,N_elements,basis_type_trial_temp,0,0,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A14 = assemble_matrix_volume_int_2D('function_mu_alpha',P,T,Tb_trial_temp,Tb_test_disp,number_of_lb_temp,number_of_lb_disp,N_basis_temp,N_basis_disp,N_elements,basis_type_trial_temp,0,0,basis_type_test_disp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A15_1 =assemble_matrix_volume_int_2D('function_K',P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,N_basis_temp,N_basis_temp,N_elements,basis_type_trial_temp,1,0,basis_type_test_temp,1,0, ref_Gauss_weight,ref_Gauss_nodes);
A15_2 =assemble_matrix_volume_int_2D('function_K',P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,N_basis_temp,N_basis_temp,N_elements,basis_type_trial_temp,0,1,basis_type_test_temp,0,1, ref_Gauss_weight,ref_Gauss_nodes);
A_t = A15_1+A15_2;
% displacement matrix
A_u = [A1+2*A2+A3 A4+A5;A8+A9 A10+A11+2*A12];
% Coupling matrix
A_u_t = [-3*A6-2*A7;-3*A13-2*A14];

%% Define Mass Matrices 
C_mat_33 = assemble_matrix_volume_int_2D('function_density_spec_ht',P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,N_basis_temp,N_basis_temp,N_elements,basis_type_trial_temp,0,0,basis_type_test_temp,0,0, ref_Gauss_weight,ref_Gauss_nodes);
M_mat_u1 = assemble_matrix_volume_int_2D('function_density',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,0,basis_type_test_disp,0,0, ref_Gauss_weight,ref_Gauss_nodes);
M_mat_u2 = assemble_matrix_volume_int_2D('function_density',P,T,Tb_trial_disp,Tb_test_disp,number_of_lb_disp,number_of_lb_disp,N_basis_disp,N_basis_disp,N_elements,basis_type_trial_disp,0,0,basis_type_test_disp,0,0, ref_Gauss_weight,ref_Gauss_nodes);
damp_mat_C = C_mat_33;
mass_mat_u = blkdiag(M_mat_u1,M_mat_u2); 
%% Initial Conditions (IC) for velocity variables: _disp_1, _disp_2 and temp variable: _temp 
     IC_disp_1 = generate_initial_vector('function_initial_conditions_disp_1',Pb_disp,0); 
     IC_disp_2 = generate_initial_vector('function_initial_conditions_disp_2',Pb_disp,0);
     IC_temp = generate_initial_vector('function_initial_conditions_temp',Pb_temp,0);
     xm_u = [IC_disp_1;IC_disp_2];
     xm_t = IC_temp;
%% Disp Solution at dt
    u1_disp_1 = generate_step_1_vector('function_initial_conditions_disp_1',Pb_disp,dt);
    u2_disp_2 = generate_step_1_vector('function_initial_conditions_disp_2',Pb_disp,dt);
    xm_plus_one_u = [u1_disp_1;u2_disp_2];


%Deal with Robin boundary condition for matrix. 
%If we don't have Robin boundary condition at all, we can comment the following line to save time and memory.
%A=treat_Robin_boundary_for_matrix_triangle_inner('function_convec_coeff',A,boundary_edges_heat,P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,N_basis_disp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_trial_temp,0,0,basis_type_test_temp,0,0);
A_t=treat_Robin_boundary_for_matrix_triangle_outer('function_convec_coeff',A_t,boundary_edges_heat,P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_trial_temp,0,0,basis_type_test_temp,0,0);
A_t=treat_Robin_boundary_for_matrix_triangle_inner('function_convec_coeff',A_t,boundary_edges_heat,P,T,Tb_trial_temp,Tb_test_temp,number_of_lb_temp,number_of_lb_temp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_trial_temp,0,0,basis_type_test_temp,0,0);


               
%% Finite Difference Scheme  
    fixed_mat = -mass_mat_u/dt^2;
    sub_fixed_mat = 2*mass_mat_u/dt^2; %For Implicit Finite Difference Scheme
    A_tilde_u = mass_mat_u/dt^2 + A_u;   %For Implicit Finite Difference Scheme
%     
%     fixed_mat = -(mass_mat_u/dt^2 + A_u/4);
%     sub_fixed_mat = 2*(mass_mat_u/dt^2 - A_u/4); 
%     A_tilde_u = (mass_mat_u/dt^2 + A_u/4); 
    A_theta_tilde_t = damp_mat_C/dt + theta*A_t;
    fixed_mass_t = damp_mat_C/dt - (1-theta)*A_t;

%% Begin Iterations
M = (end_time - start_time)/dt;
for m=0:M-1
    tm = m*dt;
    tm_current_t = (m+1)*dt;
    tm_plus_one_u = (m+1)*dt;
    tm_current_u = (m+2)*dt;
  %% Obtaining Temperatures  
                b5 = assemble_vector_b_test_2D_time('function_H',P,T,Tb_test_temp,number_of_lb_temp,tm,basis_type_test_temp,0,0,N_basis_temp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
                b_tm = b5;
                b55 = assemble_vector_b_test_2D_time('function_H',P,T,Tb_test_temp,number_of_lb_temp,tm_current_t,basis_type_test_temp,0,0,N_basis_temp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
                b_tm_plus_one = b55;
                b = theta*b_tm_plus_one + (1-theta)*b_tm;
                b_theta_tilde_t = b + fixed_mass_t*xm_t;
 % Treat Nueumann BCs
    b_theta_tilde_t = treat_Robin_boundary_for_vector_time_triangle_inner('function_qn_tilde_Robin_inner',tm_current_t,tm,theta,b_theta_tilde_t,boundary_edges_heat,P,T,Tb_test_temp,number_of_lb_temp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_test_temp,0,0);
    b_theta_tilde_t = treat_Robin_boundary_for_vector_time_triangle_outer('function_qn_tilde_Robin_outer',tm_current_t,tm,theta,b_theta_tilde_t,boundary_edges_heat,P,T,Tb_test_temp,number_of_lb_temp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_test_temp,0,0);

       % Treat Dirichlet BCs
      [A_theta_tilde_t,b_theta_tilde_t] = Treatment_of_Dirichlet_BCs_2D_time_temp_for_scheme_2('rim_temperature',A_theta_tilde_t,b_theta_tilde_t,Pb_temp,boundary_nodes_basis_heat,tm_current_t);
      % Obtain Temperature Solution 
       xm_plus_1_t = A_theta_tilde_t\b_theta_tilde_t;
 %% Obtaining Displacements                       
         b1 = assemble_vector_b_test_2D_time('function_f1',P,T,Tb_test_disp,number_of_lb_disp,tm_plus_one_u,basis_type_test_disp,0,0,N_basis_disp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
         b2 = assemble_vector_b_test_2D_time('function_b2',P,T,Tb_test_disp,number_of_lb_disp,tm_plus_one_u,basis_type_test_disp,1,0,N_basis_disp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
         b3 = assemble_vector_b_test_2D_time('function_f2',P,T,Tb_test_disp,number_of_lb_disp,tm_plus_one_u,basis_type_test_disp,0,0,N_basis_disp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
         b4 = assemble_vector_b_test_2D_time('function_b2',P,T,Tb_test_disp,number_of_lb_disp,tm_plus_one_u,basis_type_test_disp,0,1,N_basis_disp,ref_Gauss_weight,ref_Gauss_nodes,N_elements);
         b = [(b1-b2);(b3-b4)];
         b_tilde_u = b + sub_fixed_mat*xm_plus_one_u + fixed_mat*xm_u - A_u_t*xm_plus_1_t;
    
 %% Specify Boundary Type. Note that to reduce computational time one can comment out a boundary type not needed for the problem being solved

% Treat Stress BC
%       b_tilde_u = treat_stress_boundary_time_triangle_inner('function_qn_tilde_inner','function_qt_tilde_inner',tm_plus_one_u,tm,theta,b_tilde_u,boundary_edges_stress,P,T,Tb_test_disp,number_of_lb_disp,N_basis_disp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_test_disp,0,0);
      b_tilde_u = treat_stress_boundary_time_triangle_outer('function_qn_tilde_outer','function_qt_tilde_outer',tm_plus_one_u,tm,theta,b_tilde_u,boundary_edges_stress,P,T,Tb_test_disp,number_of_lb_disp,N_basis_disp,ref_Gauss_weight_1D,ref_Gauss_nodes_1D,basis_type_test_disp,0,0);

% Treat Dirichlet BCs
    [A_tilde_u,b_tilde_u] = Treatment_of_Dirichlet_BCs_2D_time_for_scheme_2('function_g1','function_g2',A_tilde_u,b_tilde_u,Pb_disp,boundary_nodes_basis_disp,tm_current_u);
% Obtain Displacement Solution Vector
    xm_plus_two_u = A_tilde_u\b_tilde_u;
% Update Vectors 
    xm_t = xm_plus_1_t;           
    xm_u = xm_plus_one_u;
    xm_plus_one_u = xm_plus_two_u;           
    
end

     u1h =   xm_plus_one_u(1:N_basis_disp,:);
     u2h =   xm_plus_one_u(N_basis_disp+1:2*N_basis_disp,:);
     uh = sqrt(u1h.^2 + u2h.^2); % u magnitude
     temp_h = xm_t;
 %Output the datum for tecplot
% 'w' specifies that it will be written. 
% Similarly 'r' is for reading and 'a' is for appending
    f1=fopen('P245_R17_tire.dat','w'); 
    fprintf(f1,'TITLE = ""\n'); 
    fprintf(f1,'VARIABLES = "X" "Y" "U, Magnitude" "U1" "U2" "Theta" \n'); 
    total_number_of_nodes=length(Pb_disp);
    total_number_of_elements=size(Tb,2);
    fprintf(f1,'ZONE N=%d,E=%d\n',total_number_of_nodes,total_number_of_elements); 
    fprintf(f1,'DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n');
for i=1:total_number_of_nodes
   fprintf(f1,'%f\t%f\t%f\t%f\t%f\t%f\n',Pb_disp(1,i),Pb_disp(2,i),uh(i),u1h(i),u2h(i),temp_h(i));  
end
fprintf('\n');
for n=1:total_number_of_elements
   fprintf(f1,'%d\t%d\t%d\n',Tb_test_disp(1,n),Tb_test_disp(2,n),Tb_test_disp(3,n)); 
end
fclose(f1);
end
