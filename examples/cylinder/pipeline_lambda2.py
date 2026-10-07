# script-version: 2.0
# Catalyst state generated using paraview version 6.0.1
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 0

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
renderView1.Set(
    ViewSize=[3040, 1558],
    CenterOfRotation=[10.000000000000004, 0.0, 3.000000000000001],
    CameraPosition=[-2.291749575748968, 4.200480496715104, 25.08686719789688],
    CameraFocalPoint=[9.546272970188394, -0.02751341719617567, 2.764589247573702],
    CameraViewUp=[0.04751987360769358, 0.9857275146611646, -0.16150271964345544],
    CameraFocalDisk=1.0,
    CameraParallelScale=25.179356624028348,
    UseColorPaletteForBackground=0,
    Background=[1.0, 1.0, 1.0],
    OSPRayMaterialLibrary=materialLibrary1,
)

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VTKHDF Reader'
fieldvtkhdf = VTKHDFReader(registrationName='field.vtkhdf', FileName=['/scratch/baconnet/externals/pySEMTools-plugins/examples/insitu_streaming/field.vtkhdf'])
fieldvtkhdf.PointArrayStatus = ['lambda2', 'mag', 'v']

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=fieldvtkhdf)
contour2.Set(
    ContourBy=['POINTS', 'lambda2'],
    ComputeNormals=0,
    GenerateTriangles=0,
    Isosurfaces=[-1.0, -0.775, -0.55, -0.32499999999999996, -0.15],
)

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=fieldvtkhdf)
contour1.Set(
    ContourBy=['POINTS', 'mag'],
    ComputeNormals=0,
    GenerateTriangles=0,
    Isosurfaces=[1e-09],
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Set(
    Representation='Surface',
    AmbientColor=[0.8313725590705872, 0.8313725590705872, 0.8313725590705872],
    ColorArrayName=['POINTS', ''],
    DiffuseColor=[0.8313725590705872, 0.8313725590705872, 0.8313725590705872],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [9.999999717180685e-10, 0.0, 0.5, 0.0, 1.0002273453935118e-09, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [9.999999717180685e-10, 0.0, 0.5, 0.0, 1.0002273453935118e-09, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'v'
vLUT = GetColorTransferFunction('v')
vLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=-0.8,
        range_max=0.864607572555542,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
contour2Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'v'],
    LookupTable=vLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-7.917351722717285, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-7.917351722717285, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vLUT in view renderView1
vLUTColorBar = GetScalarBar(vLUT, renderView1)
vLUTColorBar.Set(
    AutoOrient=0,
    Orientation='Horizontal',
    WindowLocation='Any Location',
    Position=[0.2847836538461538, 0.16600770218228506],
    Title='v',
    ComponentTitle='',
    TitleColor=[0.0, 0.0, 0.0],
    LabelColor=[0.0, 0.0, 0.0],
    ScalarBarLength=0.3300000000000004,
)

# set color bar visibility
vLUTColorBar.Visibility = 1

# show color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'v'
vPWF = GetOpacityTransferFunction('v')
vPWF.Set(
    Points=[-0.8, 0.0, 0.5, 0.0, 0.864607572555542, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.Set(
    FileName='RenderView1_{timestep:06d}{camera}.png',
    ImageResolution=[6080, 3116],
    Format='PNG',
)

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = '/scratch/baconnet/software/neko/examples/cylinder/results_insitu/catalyst_results'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
