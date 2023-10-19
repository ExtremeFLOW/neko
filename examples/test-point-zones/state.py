# state file generated using paraview version 5.11.1
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [2977, 1848]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-2.4216501447746586, 1.9053823306552697, 1.3276553382458811]
renderView1.CameraFocalPoint = [0.4999999999999997, 0.49999999999999983, 0.5]
renderView1.CameraViewUp = [0.40634303236271485, 0.9074883974357576, -0.10653707603343666]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [0.7607843137254902, 0.7607843137254902, 0.7607843137254902]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(2977, 1848)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItNek5000Reader'
field0nek5000 = VisItNek5000Reader(registrationName='field0.nek5000', FileName='/scratch/baconnet/software/neko/examples/test-point-zones/field0.nek5000')
field0nek5000.Meshes = ['mesh']
field0nek5000.PointArrays = ['temperature']

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=field0nek5000)
clip1.ClipType = 'Scalar'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'temperature']
clip1.Value = 9.9
clip1.Invert = 0

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.5, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from field0nek5000
field0nek5000Display = Show(field0nek5000, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
field0nek5000Display.Representation = 'Wireframe'
field0nek5000Display.ColorArrayName = [None, '']
field0nek5000Display.Opacity = 0.14
field0nek5000Display.SelectTCoordArray = 'None'
field0nek5000Display.SelectNormalArray = 'None'
field0nek5000Display.SelectTangentArray = 'None'
field0nek5000Display.OSPRayScaleArray = 'temperature'
field0nek5000Display.OSPRayScaleFunction = 'PiecewiseFunction'
field0nek5000Display.SelectOrientationVectors = 'None'
field0nek5000Display.ScaleFactor = 0.1
field0nek5000Display.SelectScaleArray = 'None'
field0nek5000Display.GlyphType = 'Arrow'
field0nek5000Display.GlyphTableIndexArray = 'None'
field0nek5000Display.GaussianRadius = 0.005
field0nek5000Display.SetScaleArray = ['POINTS', 'temperature']
field0nek5000Display.ScaleTransferFunction = 'PiecewiseFunction'
field0nek5000Display.OpacityArray = ['POINTS', 'temperature']
field0nek5000Display.OpacityTransferFunction = 'PiecewiseFunction'
field0nek5000Display.DataAxesGrid = 'GridAxesRepresentation'
field0nek5000Display.PolarAxes = 'PolarAxesRepresentation'
field0nek5000Display.ScalarOpacityUnitDistance = 0.04948716593053934
field0nek5000Display.OpacityArrayName = ['POINTS', 'temperature']
field0nek5000Display.SelectInputVectors = [None, '']
field0nek5000Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
field0nek5000Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
field0nek5000Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]

# show data from clip1
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# get 2D transfer function for 'temperature'
temperatureTF2D = GetTransferFunction2D('temperature')

# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.TransferFunction2D = temperatureTF2D
temperatureLUT.RGBPoints = [9.699999809265137, 0.231373, 0.298039, 0.752941, 9.849999904632568, 0.865003, 0.865003, 0.865003, 10.0, 0.705882, 0.0156863, 0.14902]
temperatureLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'temperature'
temperaturePWF = GetOpacityTransferFunction('temperature')
temperaturePWF.Points = [9.699999809265137, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]
temperaturePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'temperature']
clip1Display.LookupTable = temperatureLUT
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'temperature'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.06003848314285279
clip1Display.SelectScaleArray = 'temperature'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'temperature'
clip1Display.GaussianRadius = 0.003001924157142639
clip1Display.SetScaleArray = ['POINTS', 'temperature']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'temperature']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = temperaturePWF
clip1Display.ScalarOpacityUnitDistance = 0.06917648482996253
clip1Display.OpacityArrayName = ['POINTS', 'temperature']
clip1Display.SelectInputVectors = [None, '']
clip1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [9.699999809265137, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [9.699999809265137, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for temperatureLUT in view renderView1
temperatureLUTColorBar = GetScalarBar(temperatureLUT, renderView1)
temperatureLUTColorBar.Title = 'temperature'
temperatureLUTColorBar.ComponentTitle = ''

# set color bar visibility
temperatureLUTColorBar.Visibility = 1

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(field0nek5000)
# ----------------------------------------------------------------


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')