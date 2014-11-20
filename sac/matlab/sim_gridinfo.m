classdef sim_gridinfo
    %Class to process and store grid information used and generated by sac/smaug MHD codes
    %   Detailed explanation goes here
    
    properties
        ndimensions=3;
        grid_dimensions=ones(3,1);
        grid_left_index=0;
        %grid_right_index=0;
        grid_level=0;
        grid_parent_id=0;
        grid_particle_count=0;
        
        
        
    end
    
    methods
        
        function newobj=read_gridinfo_h5( obj, filename)
            
            newobj=sim_gridinfo;
            
            newobj.grid_dimensions=h5read(filename,'/grid_dimensions');
            newobj.grid_left_index=h5read(filename,'/grid_left_index');
            %newobj.grid_right_index=h5read(filename,'/grid_right_index');
            newobj.grid_level=h5read(filename,'/grid_level');
            newobj.grid_parent_id=h5read(filename,'/grid_parent_id');
            newobj.grid_particle_count=h5read(filename,'/grid_particle_count');
            
        end  
        
        
        function newobj=write_gridinfo_h5( obj, filename)
            
            newobj=sim_gridinfo;
            
            newobj.grid_dimensions=h5read(filename,'/grid_dimensions');
            newobj.grid_left_index=h5read(filename,'/grid_left_index');
            %newobj.grid_right_index=h5read(filename,'/grid_right_index');
            newobj.grid_level=h5read(filename,'/grid_level');
            newobj.grid_parent_id=h5read(filename,'/grid_parent_id');
            newobj.grid_particle_count=h5read(filename,'/grid_particle_count');
            
        end 
    end
    
end

