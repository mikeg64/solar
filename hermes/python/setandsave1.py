# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 
image_plane_widget1 = engine.scenes[0].children[0].children[0].children[1]
image_plane_widget1.ipw.origin = array([  0.5,   0.5,  51. ])
image_plane_widget1.ipw.point1 = array([ 128.5,    0.5,   51. ])
image_plane_widget1.ipw.point2 = array([   0.5,  128.5,   51. ])
image_plane_widget1.ipw.enabled = False
image_plane_widget1.ipw.origin = array([  0.5,   0.5,  51. ])
image_plane_widget1.ipw.point1 = array([ 128.5,    0.5,   51. ])
image_plane_widget1.ipw.point2 = array([   0.5,  128.5,   51. ])
image_plane_widget1.ipw.interactor = None
module_manager = engine.scenes[0].children[0].children[0]
module_manager.children[1:2] = []
streamline = engine.scenes[0].children[1].children[0].children[0].children[0]
streamline.seed.widget = streamline.seed.widget_list[2]
streamline.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamline.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamline.seed.widget.enabled = False
streamline.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamline.seed.widget.handle_direction = array([ 1.,  0.,  0.])
streamline.seed.widget.interactor = None
# streamline.seed.widget = <tvtk.tvtk_classes.plane_widget.PlaneWidget object at 0x136851d70>
streamline.seed.widget.origin = array([ 64.5 ,  32.75,  32.75])
streamline.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamline.seed.widget.normal = array([ 1.,  0.,  0.])
streamline.seed.widget.point1 = array([ 64.5 ,  96.25,  32.75])
streamline.seed.widget.point2 = array([ 64.5 ,  32.75,  96.25])
streamline.seed.widget.enabled = False
streamline.seed.widget.origin = array([ 64.5 ,  32.75,  32.75])
streamline.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamline.seed.widget.normal = array([ 1.,  0.,  0.])
streamline.seed.widget.point1 = array([ 64.5 ,  96.25,  32.75])
streamline.seed.widget.point2 = array([ 64.5 ,  32.75,  96.25])
streamline.seed.widget.origin = array([ 32.75,  32.75,  64.5 ])
streamline.seed.widget.origin = array([ 32.75,  32.75,  32.75])
streamline.seed.widget.origin = array([ 32.75,  32.75,  64.5 ])
streamline.seed.widget.center = array([ 64.5,  64.5,  64.5])
streamline.seed.widget.normal = array([ 0.,  0.,  1.])
streamline.seed.widget.normal = array([ 0.,  0.,  0.])
streamline.seed.widget.normal = array([ 0.,  0.,  1.])
streamline.seed.widget.point1 = array([ 96.25,  32.75,  64.5 ])
streamline.seed.widget.point1 = array([ 96.25,  96.25,  32.75])
streamline.seed.widget.point1 = array([ 96.25,  32.75,  32.75])
streamline.seed.widget.point1 = array([ 96.25,  32.75,  64.5 ])
streamline.seed.widget.point2 = array([ 32.75,  96.25,  64.5 ])
streamline.seed.widget.point2 = array([ 32.75,  32.75,  96.25])
streamline.seed.widget.point2 = array([ 32.75,  96.25,  96.25])
streamline.seed.widget.point2 = array([ 32.75,  96.25,  64.5 ])
streamline.seed.widget.normal_to_z_axis = True
streamline.stream_tracer.start_position = array([ 0.,  0.,  0.])
streamline.stream_tracer.progress = 1.0
streamline.stream_tracer.integration_direction = 'both'
streamline.stream_tracer.start_position = array([ 0.,  0.,  0.])
streamline.stream_tracer.maximum_propagation = 500.0
# streamline.clean_filter.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x143d3c890>
streamline.tube_filter.default_normal = array([ 0.,  0.,  1.])
# streamline.tube_filter.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x143d3c890>
streamline.actor.mapper.progress = 1.0
streamline.actor.mapper.scalar_range = array([ 0.        ,  0.09823792])
# streamline.actor.mapper.input = <tvtk.tvtk_classes.poly_data.PolyData object at 0x143d3c890>
streamline.streamline_type = 'tube'
streamline.tube_filter.default_normal = array([ 0.,  0.,  1.])
streamline.tube_filter.progress = 1.0
streamline.tube_filter.vary_radius = 'vary_radius_by_scalar'
streamline.tube_filter.default_normal = array([ 0.,  0.,  1.])
streamline.tube_filter.number_of_sides = 6
streamline.tube_filter.default_normal = array([ 0.,  0.,  1.])
streamline.tube_filter.capping = True
