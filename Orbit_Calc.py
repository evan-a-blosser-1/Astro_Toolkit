"""Orbital Calculations
    This is several methods of calculating orbital 
    elements & state vectors; as well as plots. 
    (Adapted from GUI_Packages.py: Asteroids_In_MASCON1)  
    -----------------------
    Author: Evan A. Blosser                              
    Date:   Fri Sep 29 2023                     
    -----------------------
         
    Function list: (assuming import as OBC)
    ---
    
    - Classic Orbital Elements to Geocentric State Vector

    
    ```python
    >>> a0 = OBC.COE_2_Geo_State(r_p_a, 
                        R_body, 
                        mu, 
                        e, 
                        i_in,  
                        w_in,  
                        nu_in, 
                        Omega_in)
    ```
    
    - State Vector to Classical Orbital Elements
    
    ```python
    >>> Orbital_Elements = OBC.State_Vec_2_COE(a,mu) 
    ```
    
                                    
    - Basic 3D Orbit Plot
    
    ```python
    >>> OBC.Orbit_3D_Plot(a)
        plt.show()
    ```
     
     
    - Orbital Gauge Cluster
    
    ```python
    >>> Astro.Orbit_Gauge_Cluster(Orbital_Elements, 
                                    a, mu,R_body, 
                                    t_inp,Graph_time,
                                    Graph_Time_Unit)
        plt.show()
    ```  
"""
############################# Imports
#####################################
# Mathmatical!!                     
import numpy as np    
from astroquery.jplhorizons import Horizons              
# System Time                       
import time                         
from tqdm import tqdm               
# Plotting & Animation              
import matplotlib.pyplot as plt   
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker  
from matplotlib.colors import LinearSegmentedColormap                      
###################################
#%% Ephemris Lookup
def ephemeris(body, start_date, end_date, step): #Define the function and the data to use

  obj = Horizons(id=body,
                 location='500@10',
                 epochs={'start':start_date, 'stop':end_date, 'step':step})
  oe = obj.elements(refplane='ecliptic')
  return oe


#%% Rotation Matrix
def rotation_matrix(angle, axis='z'):
    """Create 3D rotation matrix"""
    c = np.cos(angle)
    s = np.sin(angle)
    if axis.lower() == 'z':
        return np.array([[c, -s, 0],
                        [s, c, 0],
                        [0, 0, 1]])
    elif axis.lower() == 'y':
        return np.array([[c, 0, s],
                        [0, 1, 0],
                        [-s, 0, c]])
    elif axis.lower() == 'x':
        return np.array([[1, 0, 0],
                        [0, c, -s],
                        [0, s, c]])
        
#%% RTN Rotation 
def rtn_rot(a, vector,choice="RTN"):
    """
    Translates a vector from RTN coordinates back to Cartesian coordinates.
    Parameters:
    a (array-like): The state vector [x, y, z, vx, vy, vz].
    rtn_vector (array-like): The vector in RTN coordinates to be translated.
    Returns:
    cartesian_vector (array-like): The translated vector in Cartesian coordinates.
    """
    r = np.array([a[0], a[1], a[2]])
    r_dot = np.array([a[3], a[4], a[5]])
    
    # Normalize the position vector to get the radial direction
    r_hat = r / np.linalg.norm(r)
    
    # Compute the transverse direction
    h = np.cross(r, r_dot)
    h_hat = h / np.linalg.norm(h)
    t_hat = np.cross(h_hat, r_hat)
    
    # Construct the rotation matrix
    rotation_matrix = np.vstack((r_hat, t_hat, h_hat)).T
    
    if choice == "CART":
        # Apply the inverse rotation matrix to the RTN vector
        cartesian_vector = rotation_matrix.T @ vector
        
        return cartesian_vector
    
    elif choice == "RTN":
        # Apply the rotation matrix to the Cartesian vector
        rtn_vector = rotation_matrix @ vector
        return rtn_vector   
        
        
        
        
#%% Simulation Progress Bar
def print_progress(current_time, total_time, bar_length=11):
    fraction = current_time / total_time
    block = int(round(bar_length * fraction))
    bar = '\u25A1' * block + '-' * (bar_length - block)
    sys.stdout.write(f'\r[{bar}] {fraction:.2%}')
    sys.stdout.flush()            
            

#%% Classic Orbital Elements to Geocentric State Vector
def COE_2_Geo_State(r_p_a, R_body, mu_body, e, i_in,  w_in,  nu_in, Omega_in):
    """
    Parameters
    ----------
    r_p_a : Float
        Altitude of Perigee of orbiting body.
    R_body : Float
        Radius of body being orbited.
    mu_body : Float
        Gravitational Parameter of body being orbited.
    e : Float
        Ecccentricity of orbiting body.
    i_in : Float
        Inclination of orbiting body.
    w_in : Float
        Argument of Perigee of orbiting body.
    nu_in : Float
        True Anomoly of orbiting body.
    Omega_in : Float
        Longitude of the ascending node of orbiting body.

    Returns
    -------
    The State Vector, "state_vector"

    """
    #######################
    # Initial Calculation #
    # Cpu Clock           #################
    Initial_Calc_Start_Time = time.time() #
    #######################################
    ######################## 
    # Angle Inputs deg2Rad #
    ###############################
    i_0  = np.deg2rad(i_in)        #
    nu_0 = np.deg2rad(nu_in)        #
    Omega_0 = np.deg2rad(Omega_in)   #
    w_0  = np.deg2rad(w_in)         #
    ################################
    #
    ######################
    # Basic Calculations #
    ############################################
    # Set altitude to radius relative to Earth #
    r_p = r_p_a + R_body                       #
    # Calculate the Semi-Latus Rectum          #
    p = r_p*(1 + e)                            #
    # Position magnitue for elements           #
    r = p / (1 + e*np.cos(nu_0))               #
    ############################################
    #
    ##############################
    # Calculate perifocal vector #
    ###########################################
    # Calculate the perifocal position vector #
    r_pf_p = r*np.cos(nu_0)                   #
    r_pf_q = r*np.sin(nu_0)                   #
    # Set perifocal position vector in array  #
    r_pf = np.array([r_pf_p, r_pf_q, 0])       #
    # Calculate perifocal velocity vector       #
    v_pf_p = np.sqrt(mu_body/p)*(-np.sin(nu_0))  #
    v_pf_q = np.sqrt(mu_body/p)*(e+np.cos(nu_0)) #
    # Set perifocal velocity vector arry        #
    v_pf = np.array([v_pf_p, v_pf_q, 0])       #
    ###########################################
    #
    ##################
    # Rotaton matrix #
    ###############################################################################
    R_11 =  np.cos(Omega_0)*np.cos(w_0) - np.sin(Omega_0)*np.sin(w_0)*np.cos(i_0) #
    R_12 = -np.cos(Omega_0)*np.sin(w_0) - np.sin(Omega_0)*np.cos(w_0)*np.cos(i_0) #
    R_13 =  np.sin(Omega_0)*np.sin(i_0)                                           #
    R_21 =  np.sin(Omega_0)*np.cos(w_0) + np.cos(Omega_0)*np.sin(w_0)*np.cos(i_0) #
    R_22 = -np.sin(Omega_0)*np.sin(w_0) + np.cos(Omega_0)*np.cos(w_0)*np.cos(i_0) #
    R_23 = -np.cos(Omega_0)*np.sin(i_0)                                          #
    R_31 =  np.sin(w_0)*np.sin(i_0)                                             #
    R_32 =  np.cos(w_0)*np.sin(i_0)                                           #
    R_33 =  np.cos(i_0)                                                     #
    # Assign Matrix ########################################################
    R_matrix = np.array([[R_11, R_12, R_13], [R_21, R_22, R_23],   # 
                         [R_31, R_32, R_33]])                      #
    ################################################################
    #
    ########################
    # Calcualte Geocentric #
    ######################################
    # Geocentric Position Vector         #
    r_geo = np.matmul(R_matrix,r_pf)     #
    # Geocentric Velocity Vector         #
    v_geo = np.matmul(R_matrix,v_pf)     #
    ######################################
    #
    #############
    # Cpu Clock #########################
    Initial_Calc_End_Time = time.time() #
    #######################################
    # Initial Calculations Execution Time ####################################
    Initial_Calc_Execution = Initial_Calc_End_Time - Initial_Calc_Start_Time #
    ##########################################################################
    Vector_message = f"""
{'-'*42}
|{'-'*13} State Vector {'-'*13}|
{'-'*42}
|   Calculations completed
|   in: {Initial_Calc_Execution:.3f} Seconds
{'-'*42}
| Position vector:  
|       {r_geo[0]:.3f} km (i)  
|       {r_geo[1]:.3f} km (j) 
|       {r_geo[2]:.3f} km (k)
{'-'*42}
| Velocity vector:
|       {v_geo[0]:.3f} km (i)  
|       {v_geo[1]:.3f} km (j) 
|       {v_geo[2]:.3f} km (k)
{'-'*42}
    """
    print(Vector_message)
    state_vector = np.array([r_geo[0],
                             r_geo[1],
                             r_geo[2],
                             v_geo[0],
                             v_geo[1],
                             v_geo[2]])
    return state_vector


#%% State Vector to Classical Orbital Elements 
def State_Vec_2_COE(a,mu):
    """Creates an array of Classical Orbital Elements
        from a state vector of position & velocity

    Args:
        a (array): state_vector
        mu_body (float): Gravitational parameter of orbited body

    Returns:
        Array: Angular Momentum, Semi-Major Axis, Eccentricity, 
                Inclination, Longitude of Ascending Node,
                Argumaent of Periapsis, & True Anomaly
    """
    ###################################
    # Solve for elements at each time #
    ###################################
    # Setting empty matrices ##############################################
    rcl, vcl = [],[]                                                      #
    r_m,v_m = [],[]                                                       #
    Omega_quad,omega_quad = [],[]                                         #
    pl,n_vec,n_mag,i,Omega,omega,nu =  [],[],[],[],[],[],[]               #
    i_calc,Omega_calc,omega_calc,nu_calc,n_calc = [],[],[],[],[]          #
    hl,Axl,el,i_degl,Omega_degl,omega_degl,nu_degl = [],[],[],[],[],[],[] #
    #######################################################################
    # - for stand alone function
    #points = np.shape(a)[0]
    #####################
    # Progress Bar
    pbar1 = tqdm(total=len(a),
            bar_format='Calculating (COE) for Plot... {l_bar}{bar:10} | {n_fmt}/{total_fmt} [Elapsed: {elapsed}  Remaining: {remaining}]', 
            colour='white'
            )
    #####################
    # Start Calculating #
    #####################
    points = len(a)
    for i in range(0,points):
        pbar1.update(1)
        rc = np.array([a[i,0],a[i,1],a[i,2]])
        rcl.append(rc)
        vc = np.array([a[i,3],a[i,4],a[i,5]])
        vcl.append(vc)
        r_m = np.linalg.norm(rc)
        v_m = np.linalg.norm(vc)
        ####################
        # Angular Momentum #
        ####################################
        # Angular momentum of 1st position #
        h = np.cross(rc,vc)                #
        h_mag = np.linalg.norm(h)          #
        ####################################
        #
        #####################
        # Semi-Latus Rectum # 
        #######################
        # 1st semi-latus rec  #
        p = h_mag**2/mu       #
        #######################
        #
        ################
        # Eccentricity #
        ############################################
        # 1st eccentricity vector                  #
        e = np.cross(vc,h)/mu - rc/r_m             #
        # 1st eccentricity vector magnitude        #
        e_mag = np.sqrt((e*e).sum())               #
        ############################################
        #
        ###################
        # Semi-Major Axis #
        #######################
        # 1st semi-major axis ###########
        Ax = h_mag**2/mu/(1 - e_mag**2) #
        #################################
        #
        #################
        # Normal Vector #
        ######################################
        # Sets a normalizing vector          #
        n_calc = np.array([0, 0, 1])         #
        # 1st normal vector                  #
        n_vec = np.cross(n_calc,h)           #
        # 1st normal vector magnitude        #
        n_mag = np.sqrt((n_vec*n_vec).sum()) #
        ######################################
        #
        ###############
        # Inclination #
        ############################
        # 1st inclination          #
        i_calc = h[2]/h_mag        #
        i      = np.arccos(i_calc) #
        i_deg  = np.rad2deg(i)     #
        ############################
        #
        ###############################
        # Longitude of Ascending Node #
        ########################################
        # 1st Longitude (Omega)                #
        Ome_calc  = n_vec[0]/n_mag             #
        Omega     = np.arccos(Ome_calc)        #
        # Quadrant check                       #
        if n_vec[1] > 0:                       #
            Omega_deg = np.rad2deg(Omega)      #
        else:                                  #
            Omega_quad = 2*np.pi - Omega       #
            Omega_deg = np.rad2deg(Omega_quad) #
        ########################################
        #
        ##################################
        # Argumaent of Periapsis (omega) #
        ############################################
        # 1st Argumaent of Periapsis (omega)       #
        ome_calc  = np.dot(n_vec,e)/(n_mag*e_mag)  #
        omega     = np.arccos(ome_calc)            #
        # Quadrant check                           #
        if e[2] > 0:                               #
            omega_deg = np.rad2deg(omega)          #
        else:                                      #
            omega_quad = 2*np.pi - omega           #
            omega_deg = np.rad2deg(omega_quad)     #
        ############################################
        #
        #####################
        # True Anomaly (nu) #
        ###########################################
        # 1st true anomaly                        #
        nu_calc = np.dot(e,rc)/(e_mag*r_m)        #
        nu      = np.arccos(nu_calc)              #
        # Quadrant check                          #
        if np.dot(rc,vc) > 0:                     #             
            nu_deg  = np.rad2deg(nu)              #
        else:                                     #
            nu_quad = 2*np.pi - nu                #
            nu_deg = np.rad2deg(nu_quad)          #
        ###########################################
        #
        # debugggg
        O_list, smallo_list = [],[]
        O_list.append(Ome_calc)
        smallo_list.append(ome_calc)
        #########
        ###########
        # Storage #
        ###########
        pl.append(p)
        hl.append(h_mag)
        Axl.append(Ax)
        el.append(e_mag)
        i_degl.append(i_deg)
        Omega_degl.append(Omega_deg)
        omega_degl.append(omega_deg)
        nu_degl.append(nu_deg)
    ################################### 
    pbar1.close() #
    ###############   
    return [hl,Axl,el,i_degl,Omega_degl,omega_degl,nu_degl]


#%% Basic 3D Orbit Plot
def Orbit_3D_Plot(a):
  """MASCON Orbital plot in 3-Dimensions
      Used to plot simulated orbit around the tetrahedron center of masses,
      contianed within the asteroid's outer mesh.
  Args:
      a (array): State Vector of Orbit.
  Returns:
      plot: 3D plot of the orbit. Don't forget to call the plot with `plt.show()` !! 
  """
  ############
  # Settings #
  #####################
  # Colors
  grid_col   = '#0200FF'
  Space      = "#000000"
  orbit_line = "#F70101"
  # Grid Color                            
  plt.rcParams['grid.color'] = grid_col   
  #####################################
  # Set plot                              
  fig = plt.figure('Orbit')               
  # Set axis                              
  axis = plt.axes(projection='3d')           
  # Set Window Size                       
  fig.tight_layout()              
  # Plot assumed center                   
  cm_x = 0
  cm_y = 0
  cm_z = 0
  axis.scatter3D(cm_x,cm_y,cm_z,
                  marker='o',
                  color='#D41159')
  ##########################################
  ############################## Set Aspect
  axis.set_box_aspect([1,1,1])
  #########################################
  # Set x data to position, i. given by a #
  xline = a[:,0]                          #
  # Set y data to position, j. given by a #
  yline = a[:,1]                          #
  # Set z data to position, k. given by a #
  zline = a[:,2]                          #
  # Plot line and asteroid                #
  axis.plot3D(xline, yline, zline,        #
          color=orbit_line )              #
  # Axis Labels                           #
  axis.set_xlabel('x (km)')                 #
  axis.set_ylabel('y (km)')                  #
  axis.set_zlabel('z (km)')                   #
  axis.tick_params(axis='x', colors=grid_col) #
  axis.tick_params(axis='y', colors=grid_col) #
  axis.tick_params(axis='z', colors=grid_col) #
  axis.yaxis.label.set_color(grid_col)        #  
  axis.xaxis.label.set_color(grid_col)       #  
  axis.zaxis.label.set_color(grid_col)      #
  # Background Color                      # 
  fig.set_facecolor(Space)                #
  axis.set_facecolor(Space)                #
  # Grid Pane Color/set to clear          #
  axis.xaxis.set_pane_color((0.0, 0.0,     #
                              0.0, 0.0))  #
  axis.yaxis.set_pane_color((0.0, 0.0,     #
                              0.0, 0.0))  #
  axis.zaxis.set_pane_color((0.0, 0.0,     #
                              0.0, 0.0))  # 
  #########################################
  return


#%% Animate 3D Comparison of Two Datasets
def ani_3D_Comp(data1, data2,OBJ_File,gamma, labels=['Dataset 1', 'Dataset 2'],
                    Mesh_color='black',M_line=0.5,M_alpha=0.05,F_T=1000):
    """_summary_

    Args:
        data1 (_type_): _description_
        data2 (_type_): _description_
        labels (list, optional): _description_. Defaults to ['Dataset 1', 'Dataset 2'].

    Returns:
        _type_: _description_
        
        
    Example: AP.animate_3d_data(pos_data[:, 3:], target_data[:, 3:], ['Actual', 'Target'])
    """
    ##############################################################
    # Setup the figure and 3D axis
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    # Initialize lines with dotted style and markers
    line1, = ax.plot([], [], [], 'b:', label=labels[0], markevery=[-1], marker='o', markersize=8)
    line2, = ax.plot([], [], [], 'r:', label=labels[1], markevery=[-1], marker='o', markersize=8)
    # Calculate bounds with 20% padding
    max_range = np.max([
        np.max(data1[0]) - np.min(data1[0]),
        np.max(data1[1]) - np.min(data1[1]),
        np.max(data1[2]) - np.min(data1[2])
    ])

    # Find center point
    mid_x = (np.max(data1[0]) + np.min(data1[0])) * 0.5
    mid_y = (np.max(data1[1]) + np.min(data1[1])) * 0.5
    mid_z = (np.max(data1[2]) + np.min(data1[2])) * 0.5
    
    # Set equal limits with padding
    padding = max_range * 0.6
    ax.set_xlim(mid_x - padding, mid_x + padding)
    ax.set_ylim(mid_y - padding, mid_y + padding)
    ax.set_zlim(mid_z - padding, mid_z + padding)
    
    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal', 'box')
    ###################
    # Space Theme it for the kiddos 
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black') 
    # Add frame skipping - only animate every nth frame
    frame_skip = 2  # Adjust this value to skip more frames
    
    # Create frames array with skipped indices
    frames = range(0, len(data1[0]), frame_skip)

    def update(frame):
        # Update data for both lines
        line1.set_data(data1[0][:frame], data1[1][:frame])
        line1.set_3d_properties(data1[2][:frame])
        
        line2.set_data(data2[0][:frame], data2[1][:frame])
        line2.set_3d_properties(data2[2][:frame])
        # Ensure markers stay at end points
        line1.set_markevery([-1])  # Show marker only at latest point
        line2.set_markevery([-1])
            
        # ax.view_init(elev=30, azim=frame % 360)
        return line1, line2
    
    # Create animation
    anim = animation.FuncAnimation(
        fig, update, frames=frames, 
        interval=F_T,
        blit=True,               # Enable blitting
        cache_frame_data=False  # Disable frame caching
    )
    
    
    plt.legend()
    plt.show()


#%% Animate 3D Trajectory Data
def ani_3D_data(data1, labels=['Orbit'], F_T=0.1):
    """Animate 3D trajectory data with asteroid mesh.

    Args:
        data1: Trajectory data array (3xN)
        OBJ_File: Path to asteroid obj file
        gamma: Scaling factor for asteroid mesh
        labels: Label for trajectory. Defaults to ['Dataset 1']
        Mesh_color: Color of mesh edges. Defaults to 'black'
        M_line: Mesh line width. Defaults to 0.5
        M_alpha: Mesh transparency. Defaults to 0.05
        F_T: Animation frame time in ms. Defaults to 1000

    Example: 
        AP.ani_3D_data(pos_data[:, 3:], 'asteroid.obj', 1.0, ['Trajectory'])
    """

    # Setup figure and axis
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    line1, = ax.plot([], [], [], 'b:', label=labels[0], markevery=[-1], marker='o', markersize=2)

    # Calculate view bounds
    max_range = np.max([
        np.max(data1[0]) - np.min(data1[0]),
        np.max(data1[1]) - np.min(data1[1]),
        np.max(data1[2]) - np.min(data1[2])
    ])
    mid_x = (np.max(data1[0]) + np.min(data1[0])) * 0.5
    mid_y = (np.max(data1[1]) + np.min(data1[1])) * 0.5
    mid_z = (np.max(data1[2]) + np.min(data1[2])) * 0.5
    
    # Set view limits
    padding = max_range * 0.6
    ax.set_xlim(mid_x - padding, mid_x + padding)
    ax.set_ylim(mid_y - padding, mid_y + padding)
    ax.set_zlim(mid_z - padding, mid_z + padding)
    
    # Style settings
    ax.set_aspect('equal', 'box')
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

    # Animation settings
    frame_skip = 2
    frames = range(0, len(data1[0]), frame_skip)

    def update(frame):
        line1.set_data(data1[0][:frame], data1[1][:frame])
        line1.set_3d_properties(data1[2][:frame])
        line1.set_markevery([-1])
        return (line1,)
    
    anim = animation.FuncAnimation(
        fig, update, frames=frames,
        interval=F_T, blit=True,
        cache_frame_data=False
    )
    
    plt.legend()
    plt.show()

#%% Animate multiple 3D trajectories
def ani_3D_NOrbs( *trajectories, labels=None,  F_T=1000):
    """Animate multiple 3D trajectories with asteroid mesh.
    Expects trajectories in Nx3 format (no transpose needed).
    """
    if labels is None:
        labels = [f'Trajectory {i+1}' for i in range(len(trajectories))]
    
    colors = plt.cm.rainbow(np.linspace(0, 1, len(trajectories)))
    
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    lines = []
    for i, color in enumerate(colors):
        line, = ax.plot([], [], [], ':', label=labels[i], 
                       markevery=[-1], marker='o', markersize=8,
                       color=color)
        lines.append(line)
    
    # Calculate bounds using first trajectory (now Nx3)
    first_traj = trajectories[0]
    max_range = np.max([
        np.max(first_traj[:, 0]) - np.min(first_traj[:, 0]),
        np.max(first_traj[:, 1]) - np.min(first_traj[:, 1]),
        np.max(first_traj[:, 2]) - np.min(first_traj[:, 2])
    ])
    mid_x = (np.max(first_traj[:, 0]) + np.min(first_traj[:, 0])) * 0.5
    mid_y = (np.max(first_traj[:, 1]) + np.min(first_traj[:, 1])) * 0.5
    mid_z = (np.max(first_traj[:, 2]) + np.min(first_traj[:, 2])) * 0.5
    
    padding = max_range * 0.6
    ax.set_xlim(mid_x - padding, mid_x + padding)
    ax.set_ylim(mid_y - padding, mid_y + padding)
    ax.set_zlim(mid_z - padding, mid_z + padding)
    
    ax.set_aspect('equal', 'box')
    ax.grid(False)
    ax.set_axis_off()
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

    frame_skip = 2
    max_frames = max(len(traj) for traj in trajectories)
    frames = range(0, max_frames, frame_skip)

    def update(frame):
        for line, traj in zip(lines, trajectories):
            if frame < len(traj):
                line.set_data(traj[:frame, 0], traj[:frame, 1])
                line.set_3d_properties(traj[:frame, 2])
                line.set_markevery([-1])
        return tuple(lines)
    
    anim = animation.FuncAnimation(
        fig, update, frames=frames,
        interval=F_T, blit=True,
        cache_frame_data=False
    )
    
    plt.legend()
    plt.show()



#%% Orbital Gauge Cluster
def Orbit_Gauge_Cluster(Orbital_Elements, a, mu,R_body, t_inp,Graph_time,Graph_Time_Unit):
  """Orbital Gauge CLuster
      This is a sub-plot of orbital elements along with plot limits.

  Args:
      Orbnital_Elements (Pandas DataFrame): A DataFrame of Orbital Elements to be plotted. 
      a (array): ODEINT output in the form of the state vector.
      mu (flaot): Gravitational Parameter of body being orbited.
      R_body (flaot): Radius of body being orbited.
      t_inp (flaot): Duration of plot.
      Graph_time (flaot): Conversion to seconds selected.
      Graph_Time_Unit (string): The unit of time for x-axis label.

  Returns:
        Sub-Plots: Orbital-Gauge-Cluster! Don't forget to call the plot with `plt.show()` !!
  """
  ##########
  # Colors #
  ##########
  # COE Output Plots
  Tick_Mark_Col    = '#1A85FF'
  major_grid_col   = "#00E2F9"
  minor_grid_col   = "#00C8DC"
  Background       = "#000000"
  Font_color       = "#02FF1F"
  Plot_line_col    = "#f2ffb7"
  Limit_Line_Color = '#ff06b5'
  ###########################
  # Fromat Graph_Time 
  n = t_inp*Graph_time
  t_span  = np.linspace(0,n,n)
  Formated_Graph_Time = t_span/Graph_time
  #############################
  # Orbital Elements Plotting #
  ################################################################
  # Set subplot                                                  
  figure1, axis = plt.subplots( 3, 3,  figsize=(12, 10),         
                        facecolor=Background)   
  # Figure MAIN Title
  figure1.suptitle(f'Orbit Simulated: {t_inp} {Graph_Time_Unit}',
                  fontsize=16, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                
  # Window Title                                                 
  figure1.canvas.manager.set_window_title(                       
      'Orbital Gauge Cluster')                                 
  #####################
  # Angular Momentum  ################
  axis[0, 0].plot(Formated_Graph_Time,Orbital_Elements[0], color=Plot_line_col) 
  # Tick Settings  
  axis[0, 0].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[0, 0].yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1)) 
  # Labels & Colors         
  axis[0, 0].set_title("Angular Momentum (h)",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[0, 0].set_xlabel('Time ({})'.format(Graph_Time_Unit) ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[0, 0].set_ylabel('\u0394h (km^2/s)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[0, 0].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[0, 0].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[0, 0].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[0, 0].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)
  axis[0, 0].set_facecolor(Background) 
  ###############################
  #               
  ###########################
  # Semi-Major Axis Plot    ###############
  axis[0, 1].plot(Formated_Graph_Time,Orbital_Elements[1], color=Plot_line_col)
  # Tick Settings  
  axis[0, 1].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[0, 1].yaxis.set_minor_locator(ticker.MultipleLocator(base=100))     
  # Labels & Colors            
  axis[0, 1].set_title("Semi-Major Axis (A)",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[0, 1].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[0, 1].set_ylabel('\u0394A (km)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[0, 1].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[0, 1].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[0, 1].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[0, 1].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)                  
  axis[0, 1].set_facecolor(Background)  
  #################################
  #
  ######################
  # Eccentricity Plot  ################
  #####################################
  # Set limit at e = 0.05
  Eccentricity_Limit = 0.05* np.ones_like(Formated_Graph_Time+30)
  # Plot Limit
  axis[0, 2].plot(Formated_Graph_Time, Eccentricity_Limit,linestyle='dashdot',
                  color=Limit_Line_Color,label='Ecc. Lim')
  # Set Legend
  Leg_Alt = axis[0, 2].legend(loc='upper right', bbox_to_anchor=(1.2, 1.15))
  Leg_Alt.get_frame().set_facecolor(Background)  
  Leg_Alt.get_texts()[0].set_color(Font_color)  
  # Plot Data
  axis[0, 2].plot(Formated_Graph_Time,Orbital_Elements[2], color=Plot_line_col) 
  # Tick Settings  
  axis[0, 2].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[0, 2].yaxis.set_minor_locator(ticker.MultipleLocator(base=0.1))  
  # Labels & Colors             
  axis[0, 2].set_title("Eccentricity:",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[0, 2].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[0, 2].set_ylabel('\u0394e',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[0, 2].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[0, 2].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[0, 2].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[0, 2].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)                
  axis[0, 2].set_facecolor(Background)     
  ###############################
  #
  ######################
  # Inclination Plot   ######################### 
  ###############################################
  # Set limit at i = 0.001
  Inclination_Limit = 0.001 * np.ones_like(Formated_Graph_Time)
  # Plot Limit
  axis[1, 0].plot(Formated_Graph_Time, Inclination_Limit,linestyle='dashdot',
                  color=Limit_Line_Color,label='Inc. lim')
  # Set Legend
  Leg_Alt = axis[1, 0].legend(loc='upper right', bbox_to_anchor=(1.2, 1.15))
  Leg_Alt.get_frame().set_facecolor(Background)  
  Leg_Alt.get_texts()[0].set_color(Font_color)  
  # Plot Data                                  
  axis[1, 0].plot(Formated_Graph_Time,Orbital_Elements[3], color=Plot_line_col)     
  # Tick Settings  
  axis[1, 0].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[1, 0].yaxis.set_minor_locator(ticker.MultipleLocator(base=5))   
  # Labels & Colors       
  axis[1, 0].set_title("Inclination:",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[1, 0].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[1, 0].set_ylabel('\u0394i (degrees)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[1, 0].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[1, 0].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[1, 0].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[1, 0].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)                  
  axis[1, 0].set_facecolor(Background)     
  ##################################
  #
  ####################################
  # Longitude of Ascending Node Plot ##############
  axis[1, 1].plot(Formated_Graph_Time,Orbital_Elements[4], color=Plot_line_col)
  # Tick Settings  
  axis[1, 1].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[1, 1].yaxis.set_minor_locator(ticker.MultipleLocator(base=25))   
  # Labels & Colors       
  axis[1, 1].set_title("Longitude of Ascending Node (\u03A9)",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[1, 1].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[1, 1].set_ylabel('\u0394\u03A9 (degrees)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[1, 1].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[1, 1].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[1, 1].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[1, 1].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)
  axis[1, 1].set_facecolor(Background)  
  ##################################
  #  
  ########################
  # Orbital Details List #
  #############################################################
  # Period Calc ###############################################
  ax_avg            = np.mean(Orbital_Elements[1])
  ax_uncert         = (np.std(Orbital_Elements[1])/ax_avg)*100
  Orbit_Period_sec  = 2*np.pi*np.sqrt((ax_avg**3)/mu)
  Total_Sim_Time    = t_inp*Graph_time
  Number_Of_Orbits  = Total_Sim_Time/Orbit_Period_sec
  Orbit_Period_unit = 'Sec.'
  # Change Units of Period
  if Orbit_Period_sec > 120:
      Orbit_Period_min = Orbit_Period_sec/60
      Orbit_Period     = Orbit_Period_min
      Orbit_Period_unit = 'Min.'
      if Orbit_Period_min > 120:
          Orbit_Period_hr = Orbit_Period_min/60
          Orbit_Period     = Orbit_Period_hr
          Orbit_Period_unit = 'Hrs.'
          if Orbit_Period_hr > 48:
              Orbit_Period_day = Orbit_Period_hr/24
              Orbit_Period     = Orbit_Period_day
              Orbit_Period_unit = 'Days' 
  ############################################################
  # Orbit Energy #############################################
  Orb_ecc       = np.mean(Orbital_Elements[2])
  ecc_uncert    = (np.std(Orbital_Elements[2])/Orb_ecc)*100
  Energy_uncert = ax_uncert + ecc_uncert
  Orb_apogee    = ax_avg*(1 + Orb_ecc)
  Orb_perigee   = ax_avg*(1 - Orb_ecc)
  Orbit_Energy  = - mu/(Orb_apogee + Orb_perigee)
  ############################################################
  # Pro/Retro    #############################################
  inclin_avg = np.mean(Orbital_Elements[3])
  if inclin_avg < 90:
      Orbit_Type_Output = 'PROGRADE'  
  elif inclin_avg > 90:
      Orbit_Type_Output = 'RETROGRADE'  
  elif inclin_avg == 90:
      Orbit_Type_Output = 'POLAR'
  #######################################################################
  ########################## Text Output ################################
  Display_Text_Out = [f"""
  {'-'*42}
  {Orbit_Type_Output} orbit 
  {'-'*42}
  {'-'*42}
  Period = {Orbit_Period:.3f} {Orbit_Period_unit} 
  (+/-) {ax_uncert:.2f} %
  {'-'*42}
  Number of simulated orbits: {Number_Of_Orbits:.2f} 
  {'-'*42}
  Energy = {Orbit_Energy:.3e} km^2/s^2 
  (+/-) {Energy_uncert:.2f} % 
  {'-'*42}
  """
  ]
  axis[1, 2].text(0.5, 0.5, Display_Text_Out[0],
                              horizontalalignment='center',
                              verticalalignment='center',
                              fontsize=12, color= Font_color)
  axis[1, 2].set_facecolor(Background)
  ##################################################################
  ##################################################################

  #############################
  # Argument of Perigee Plot  ###############
  axis[2, 0].plot(Formated_Graph_Time,Orbital_Elements[5], color=Plot_line_col)    
  # Tick Settings  
  axis[2, 0].xaxis.set_minor_locator(ticker.MultipleLocator(base=.5)) 
  axis[2, 0].yaxis.set_minor_locator(ticker.MultipleLocator(base=25))    
  # Labels & Colors        
  axis[2, 0].set_title("Argument of Perigee (\u03C9)",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[2, 0].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[2, 0].set_ylabel('\u0394\u03C9 (degrees)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[2, 0].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[2, 0].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[2, 0].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[2, 0].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)      
  axis[2, 0].set_facecolor(Background) 
  #################################
  #
  #####################
  # True Anomaly Plot ##########
  axis[2, 1].plot(Formated_Graph_Time,Orbital_Elements[6], color=Plot_line_col)     
  # Tick Settings  
  axis[2, 1].xaxis.set_minor_locator(ticker.MultipleLocator(base=0.5)) 
  axis[2, 1].yaxis.set_minor_locator(ticker.MultipleLocator(base=25))         
  # Labels & Colors            
  axis[2, 1].set_title("True Anomaly (\u03BD)",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[2, 1].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,fontsize=12, 
                        fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[2, 1].set_ylabel('\u0394\u03BD (degrees)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[2, 1].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[2, 1].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[2, 1].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[2, 1].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)
  axis[2, 1].set_facecolor(Background) 
  ####################################
  #
  ##################################
  # Altimeter!! Big Brain Idea lol ###########################
  ############################################################
  # Set limit at mean radius of asteroid, make a list of 
  # these points to graph
  Altimeter_Limit = R_body* np.ones_like(Formated_Graph_Time)
  # Plot Limit
  axis[2, 2].plot(Formated_Graph_Time, Altimeter_Limit,linestyle='dashdot',
                  color=Limit_Line_Color,label='Crash Lim')
  # Set Legend
  Leg_Alt = axis[2, 2].legend(loc='upper right', bbox_to_anchor=(1.2, 1.15))
  Leg_Alt.get_frame().set_facecolor(Background)  
  Leg_Alt.get_texts()[0].set_color(Font_color)   
  # Find Altitude:
  Altitude_Sim = (a[:,0]**2 + a[:,1]**2 + a[:,2]**2 )**(1/2)
  # Plot Data
  axis[2, 2].plot(Formated_Graph_Time,Altitude_Sim, color=Plot_line_col)     
  # Tick Settings  
  axis[2, 2].xaxis.set_minor_locator(ticker.MultipleLocator(base=1)) 
  axis[2, 2].yaxis.set_minor_locator(ticker.MultipleLocator(base=6))         
  # Labels & Colors            
  axis[2, 2].set_title("Altimeter:",
                      fontsize=14, fontweight='bold',fontdict={'family': 'Consolas'},color=Font_color)                     
  axis[2, 2].set_xlabel('Time ({})'.format(Graph_Time_Unit)  ,fontsize=12, 
                        fontweight='light',fontdict={'family': 'Consolas'},color=Font_color)                                  
  axis[2, 2].set_ylabel('Altitude (km)',
                        fontsize=12, fontweight='light',fontdict={'family': 'Consolas'},color=Font_color) 
  axis[2, 2].tick_params(axis='x', colors=Tick_Mark_Col) 
  axis[2, 2].tick_params(axis='y', colors=Tick_Mark_Col)  
  axis[2, 2].grid(which='major', linestyle='--',linewidth=0.5, color=major_grid_col)
  axis[2, 2].grid(which='minor', linestyle=':', linewidth='0.5', color=minor_grid_col)
  axis[2, 2].set_facecolor(Background) 

  ###########################################################
  #################### Combine all the operations and display                       
  plt.tight_layout()                                 
  ################################################################
  
  
  
  
  
#%% TimePrint Calcuation timer
def TimePrint(Start,End):
    Calc_Execution = End - Start
    ################################
    # Convert to Minutes if to long
    calc_t_unit = 'Seconds'
    # Minutes
    if Calc_Execution > 120:
        Calc_Execution_Min = Calc_Execution/60
        Calc_Execution = Calc_Execution_Min 
        calc_t_unit = 'Minutes'
        # Hours
        if  Calc_Execution_Min > 120:
            Calc_Execution_Hr = Calc_Execution_Min/60
            Calc_Execution = Calc_Execution_Hr
            calc_t_unit = 'Hours'
    ##############################################################
    ODE_Solved_message = f"""
    {'-'*42}
    |{'-'*9} Calculation Complete! {'-'*9}|
    {'-'*42}
    | Completed in: {Calc_Execution:.3f} {calc_t_unit}
    {'-'*42}
    """
    print(ODE_Solved_message)
    return Calc_Execution
  
  
#%% ColorBlind_Palette ()
def ColrBlind_Palette(Direction):
    """Colorblind pallete for plotting, from:
    https://davidmathlogic.com/colorblind/
    https://personal.sron.nl/~pault/
    
    Args:
        Direction (int): Reverse direction of color palette or returns 
                                the list of colors.

    Returns:
        cmap: color-blind friendly pallette :)
    """
    print(Direction)
    
    if Direction == 42:
        Color_Object = ['#882255', '#AA4499', '#CC6677', '#DDCC77', '#88CCEE', '#44AA99', '#117733', '#332288']
    
    if Direction == 0:
        Color_Palette = ['#882255', '#AA4499', '#CC6677', '#DDCC77', '#88CCEE', '#44AA99', '#117733', '#332288']
        gradient_cmap = LinearSegmentedColormap.from_list("custom_gradient",Color_Palette)
        Color_Object = gradient_cmap
    
    if Direction == 1:
        Color_P_Rev = ['#332288', '#117733', '#44AA99', '#88CCEE', '#DDCC77','#CC6677','#AA4499','#882255']
        gradient_cmap_r = LinearSegmentedColormap.from_list("custom_gradient_r",Color_P_Rev)
        Color_Object = gradient_cmap_r
        
    return Color_Object
    
#%% Find Minimum Value and Index
#################################
def Find_MinIndx(a):
    """Find the Minimum Value and Index of an Array

    Args:
        a (Numpy array): A Numpy array of any size

    Returns:
        Shape: _description_
    """
    Min_Val = np.min(a)
    print('-'*42)
    print(f'| Minimum of the input: {Min_Val}')
    print('-'*42)
    Index = np.unravel_index(a.argmin(), a.shape)
    print('-'*42)
    print(f"| The minimum's index is: {Index}")
    print('-'*42)
    return Index