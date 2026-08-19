#!/usr/bin/python3
#ReadME: This script merges the stl files from the bodiesInfo/ directory into a single file in the STLMerged/ directory.
#created by OStudenik 
import os
import numpy # as np
def canBeConvertedToFloat(input):
    try:
        float(input)
        return True
    except ValueError:
        return False

def getParticles_List(caseDir):
    "Returns a list with the particle numbers in the case directory."
    Directory = caseDir
    Strings   = list(set([numStr for numStr in os.listdir(Directory)]))
    Full_List = []
    for i in range(len(Strings)):
        string_to_save = caseDir+'/'+Strings[i]
        if(Strings[i][-4:] == '.stl'):
                Full_List.append(Strings[i][:-4])
    if(len(Full_List) > 0):            
        Full_List.sort(key = float)
    return Full_List 

Full_List = os.listdir('bodiesInfo/')
Full_List.sort(key = float)


# Implicit-sphere extension.  The original STL path below is unchanged.
# A sphere timestep is stored as one compact point/radius VTP file.  ParaView
# renders those points with one instanced Sphere source; no per-particle surface
# triangles are written by this script.
import re
import xml.etree.ElementTree as ElementTree
from xml.sax.saxutils import escape

_FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_SPHERE_BLOCK = re.compile(r"\bsphere\s*\{(.*?)\}", re.DOTALL)
_BODY_ID = re.compile(r"\bbodyId\s+(-?\d+)\s*;")
_POSITION = re.compile(r"\bposition\s*\(([^)]*)\)\s*;")
_RADIUS = re.compile(r"\bradius\s+(" + _FLOAT + r")\s*;")


def _recorded_stl_exists():
    for time_name in Full_List:
        stl_path = 'bodiesInfo/' + time_name + '/stlFiles'
        if os.path.isdir(stl_path) and len(getParticles_List(stl_path)) > 0:
            return True
    return False


def _read_sphere(info_path):
    with open(info_path, 'r') as file:
        text = file.read()

    sphere_match = _SPHERE_BLOCK.search(text)
    if sphere_match is None:
        return None

    sphere_text = sphere_match.group(1)
    body_id_match = _BODY_ID.search(text)
    position_match = _POSITION.search(sphere_text)
    radius_match = _RADIUS.search(sphere_text)

    if body_id_match is None or position_match is None or radius_match is None:
        raise ValueError('incomplete sphere record in ' + info_path)

    position = tuple(float(value) for value in position_match.group(1).split())
    if len(position) != 3:
        raise ValueError(
            'sphere position must contain three components in ' + info_path
        )

    radius = float(radius_match.group(1))
    if radius <= 0.0:
        raise ValueError('non-positive sphere radius in ' + info_path)

    return int(body_id_match.group(1)), position, radius


def _include_sphere(bounds, sphere):
    position = sphere[1]
    radius = sphere[2]
    lower = [position[direction] - radius for direction in range(3)]
    upper = [position[direction] + radius for direction in range(3)]

    if bounds is None:
        return lower + upper

    for direction in range(3):
        bounds[direction] = min(bounds[direction], lower[direction])
        bounds[direction + 3] = max(
            bounds[direction + 3], upper[direction]
        )
    return bounds


def _write_sphere_vtp(output_path, spheres):
    temporary_path = output_path + '.tmp'
    number_of_spheres = len(spheres)

    with open(temporary_path, 'w') as output:
        output.write('<?xml version="1.0"?>\n')
        output.write(
            '<VTKFile type="PolyData" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        output.write('  <PolyData>\n')
        output.write(
            '    <Piece NumberOfPoints="%d" NumberOfVerts="%d" '
            'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n'
            % (number_of_spheres, number_of_spheres)
        )
        output.write('      <PointData Scalars="radius">\n')
        output.write(
            '        <DataArray type="Float64" Name="radius" '
            'format="ascii">\n'
        )
        for body_id, position, radius in spheres:
            output.write('          %.16g\n' % radius)
        output.write('        </DataArray>\n')
        output.write(
            '        <DataArray type="Int64" Name="bodyId" '
            'format="ascii">\n'
        )
        for body_id, position, radius in spheres:
            output.write('          %d\n' % body_id)
        output.write('        </DataArray>\n')
        output.write('      </PointData>\n')
        output.write('      <CellData/>\n')
        output.write('      <Points>\n')
        output.write(
            '        <DataArray type="Float64" NumberOfComponents="3" '
            'format="ascii">\n'
        )
        for body_id, position, radius in spheres:
            output.write(
                '          %.16g %.16g %.16g\n'
                % (position[0], position[1], position[2])
            )
        output.write('        </DataArray>\n')
        output.write('      </Points>\n')
        output.write('      <Verts>\n')
        output.write(
            '        <DataArray type="Int64" Name="connectivity" '
            'format="ascii">\n          '
        )
        output.write(' '.join(str(index) for index in range(number_of_spheres)))
        output.write('\n        </DataArray>\n')
        output.write(
            '        <DataArray type="Int64" Name="offsets" '
            'format="ascii">\n          '
        )
        output.write(
            ' '.join(str(index + 1) for index in range(number_of_spheres))
        )
        output.write('\n        </DataArray>\n')
        output.write('      </Verts>\n')
        output.write('    </Piece>\n')
        output.write('  </PolyData>\n')
        output.write('</VTKFile>\n')

    os.replace(temporary_path, output_path)


def _read_vtp_bounds(input_path):
    root = ElementTree.parse(input_path).getroot()
    points_array = root.find('.//Points/DataArray')
    radius_array = None
    for data_array in root.findall('.//PointData/DataArray'):
        if data_array.get('Name') == 'radius':
            radius_array = data_array
            break

    if points_array is None or radius_array is None:
        raise ValueError('missing sphere position/radius data in ' + input_path)

    coordinates = [float(value) for value in (points_array.text or '').split()]
    radii = [float(value) for value in (radius_array.text or '').split()]
    if len(coordinates) != 3*len(radii):
        raise ValueError('inconsistent sphere data in ' + input_path)

    bounds = None
    for index, radius in enumerate(radii):
        position = tuple(coordinates[3*index:3*index + 3])
        bounds = _include_sphere(bounds, (index, position, radius))
    return bounds, len(radii)


def _merge_bounds(bounds, addition):
    if addition is None:
        return bounds
    if bounds is None:
        return addition
    for direction in range(3):
        bounds[direction] = min(bounds[direction], addition[direction])
        bounds[direction + 3] = max(
            bounds[direction + 3], addition[direction + 3]
        )
    return bounds


def _write_pvd(records):
    output_path = 'STLMerged/Sphere_Results.pvd'
    temporary_path = output_path + '.tmp'

    with open(temporary_path, 'w') as output:
        output.write('<?xml version="1.0"?>\n')
        output.write(
            '<VTKFile type="Collection" version="0.1" '
            'byte_order="LittleEndian">\n'
        )
        output.write('  <Collection>\n')
        for time_name, file_name in records:
            output.write(
                '    <DataSet timestep="%s" group="" part="0" file="%s"/>\n'
                % (escape(time_name), escape(file_name))
            )
        output.write('  </Collection>\n')
        output.write('</VTKFile>\n')

    os.replace(temporary_path, output_path)
    return output_path


def _write_paraview_state(pvd_path, bounds):
    if bounds is None:
        focal_point = (0.0, 0.0, 0.0)
        parallel_scale = 1.0
        camera_distance = 4.0
    else:
        focal_point = tuple(
            0.5*(bounds[direction] + bounds[direction + 3])
            for direction in range(3)
        )
        span = tuple(
            bounds[direction + 3] - bounds[direction]
            for direction in range(3)
        )
        parallel_scale = 1.1*max(
            0.5*span[1],
            0.5*span[0]/(4.0/3.0),
            0.05*max(span),
            1.0e-12
        )
        camera_distance = max(4.0*max(span), 4.0*parallel_scale, 1.0e-9)

    camera_position = (
        focal_point[0],
        focal_point[1],
        focal_point[2] + camera_distance
    )
    absolute_pvd_path = escape(
        os.path.abspath(pvd_path), {'"': '&quot;', "'": '&apos;'}
    )
    state_path = 'STLMerged/Sphere_Results.pvsm'
    temporary_path = state_path + '.tmp'

    state = '''<GenericParaViewApplication>
  <ServerManagerState version="5.13.1">
    <Proxy group="misc" type="TimeKeeper" id="100" servers="16">
      <Property name="TimeSources" id="100.TimeSources" number_of_elements="1"><Proxy value="200"/></Property>
      <Property name="Views" id="100.Views" number_of_elements="1"><Proxy value="400"/></Property>
    </Proxy>
    <Proxy group="animation" type="AnimationScene" id="110" servers="16">
      <Property name="Cues" id="110.Cues" number_of_elements="1"><Proxy value="120"/></Property>
      <Property name="PlayMode" id="110.PlayMode" number_of_elements="1"><Element index="0" value="2"/></Property>
      <Property name="TimeKeeper" id="110.TimeKeeper" number_of_elements="1"><Proxy value="100"/></Property>
      <Property name="ViewModules" id="110.ViewModules" number_of_elements="1"><Proxy value="400"/></Property>
    </Proxy>
    <Proxy group="animation" type="TimeAnimationCue" id="120" servers="16">
      <Property name="AnimatedPropertyName" id="120.AnimatedPropertyName" number_of_elements="1"><Element index="0" value="Time"/></Property>
      <Property name="AnimatedProxy" id="120.AnimatedProxy" number_of_elements="1"><Proxy value="100"/></Property>
      <Property name="Enabled" id="120.Enabled" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="UseAnimationTime" id="120.UseAnimationTime" number_of_elements="1"><Element index="0" value="1"/></Property>
    </Proxy>
    <Proxy group="sources" type="PVDReader" id="200" servers="1">
      <Property name="FileName" id="200.FileName" number_of_elements="1"><Element index="0" value="%(pvd)s"/></Property>
      <Property name="PointArrayStatus" id="200.PointArrayStatus" number_of_elements="4"><Element index="0" value="radius"/><Element index="1" value="1"/><Element index="2" value="bodyId"/><Element index="3" value="1"/></Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="300" servers="21">
      <Property name="Center" id="300.Center" number_of_elements="3"><Element index="0" value="0"/><Element index="1" value="0"/><Element index="2" value="0"/></Property>
      <Property name="PhiResolution" id="300.PhiResolution" number_of_elements="1"><Element index="0" value="32"/></Property>
      <Property name="Radius" id="300.Radius" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="ThetaResolution" id="300.ThetaResolution" number_of_elements="1"><Element index="0" value="32"/></Property>
    </Proxy>
    <Proxy group="views" type="RenderView" id="400" servers="21">
      <Property name="Representations" id="400.Representations" number_of_elements="1"><Proxy value="500"/></Property>
      <Property name="ViewSize" id="400.ViewSize" number_of_elements="2"><Element index="0" value="800"/><Element index="1" value="600"/></Property>
      <Property name="CenterOfRotation" id="400.CenterOfRotation" number_of_elements="3"><Element index="0" value="%(fx).16g"/><Element index="1" value="%(fy).16g"/><Element index="2" value="%(fz).16g"/></Property>
      <Property name="CameraFocalPoint" id="400.CameraFocalPoint" number_of_elements="3"><Element index="0" value="%(fx).16g"/><Element index="1" value="%(fy).16g"/><Element index="2" value="%(fz).16g"/></Property>
      <Property name="CameraPosition" id="400.CameraPosition" number_of_elements="3"><Element index="0" value="%(cx).16g"/><Element index="1" value="%(cy).16g"/><Element index="2" value="%(cz).16g"/></Property>
      <Property name="CameraViewUp" id="400.CameraViewUp" number_of_elements="3"><Element index="0" value="0"/><Element index="1" value="1"/><Element index="2" value="0"/></Property>
      <Property name="CameraParallelProjection" id="400.CameraParallelProjection" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="CameraParallelScale" id="400.CameraParallelScale" number_of_elements="1"><Element index="0" value="%(scale).16g"/></Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="500" servers="21">
      <Property name="Input" id="500.Input" number_of_elements="1"><Proxy value="200" output_port="0"/></Property>
      <Property name="Representation" id="500.Representation" number_of_elements="1"><Element index="0" value="3D Glyphs"/></Property>
      <Property name="Visibility" id="500.Visibility" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="GlyphType" id="500.GlyphType" number_of_elements="1"><Proxy value="300" output_port="0"/><Domain name="input_type" id="500.GlyphType.input_type"/><Domain name="proxy_list" id="500.GlyphType.proxy_list"><Proxy value="300"/></Domain></Property>
      <Property name="Masking" id="500.Masking" number_of_elements="1"><Element index="0" value="0"/></Property>
      <Property name="Orient" id="500.Orient" number_of_elements="1"><Element index="0" value="0"/></Property>
      <Property name="ScaleFactor" id="500.ScaleFactor" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="ScaleMode" id="500.ScaleMode" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="Scaling" id="500.Scaling" number_of_elements="1"><Element index="0" value="1"/></Property>
      <Property name="SelectScaleArray" id="500.SelectScaleArray" number_of_elements="1"><Element index="0" value="radius"/></Property>
      <Property name="SetScaleArray" id="500.SetScaleArray" number_of_elements="5"><Element index="0" value=""/><Element index="1" value=""/><Element index="2" value=""/><Element index="3" value="0"/><Element index="4" value="radius"/></Property>
      <Property name="ColorArrayName" id="500.ColorArrayName" number_of_elements="5"><Element index="0" value=""/><Element index="1" value=""/><Element index="2" value=""/><Element index="3" value="0"/><Element index="4" value=""/></Property>
      <Property name="DiffuseColor" id="500.DiffuseColor" number_of_elements="3"><Element index="0" value="0.8"/><Element index="1" value="0.8"/><Element index="2" value="0.8"/></Property>
    </Proxy>
    <Proxy group="misc" type="ViewLayout" id="600" servers="20"><Layout number_of_elements="1"><Item direction="0" fraction="0.5" view="400"/></Layout></Proxy>
    <ProxyCollection name="animation"><Item id="110" name="AnimationScene1"/><Item id="120" name="TimeAnimationCue1"/></ProxyCollection>
    <ProxyCollection name="layouts"><Item id="600" name="Layout #1"/></ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.500"><Item id="300" name="GlyphType" logname="GeometryRepresentation1/Glyph3DRepresentation/GlyphType/SphereSource"/></ProxyCollection>
    <ProxyCollection name="representations"><Item id="500" name="GeometryRepresentation1" logname="GeometryRepresentation1"/></ProxyCollection>
    <ProxyCollection name="sources"><Item id="200" name="Sphere_Results.pvd" logname="Sphere_Results.pvd"/></ProxyCollection>
    <ProxyCollection name="timekeeper"><Item id="100" name="TimeKeeper1"/></ProxyCollection>
    <ProxyCollection name="views"><Item id="400" name="RenderView1" logname="RenderView1"/></ProxyCollection>
    <CustomProxyDefinitions/>
    <Links/>
  </ServerManagerState>
</GenericParaViewApplication>
''' % {
        'pvd': absolute_pvd_path,
        'fx': focal_point[0],
        'fy': focal_point[1],
        'fz': focal_point[2],
        'cx': camera_position[0],
        'cy': camera_position[1],
        'cz': camera_position[2],
        'scale': parallel_scale
    }

    with open(temporary_path, 'w') as output:
        output.write(state)
    os.replace(temporary_path, state_path)


def _merge_implicit_spheres():
    if not os.path.isdir('STLMerged'):
        os.mkdir('STLMerged')

    records = []
    bounds = None
    number_of_spheres = 0

    for output_index, time_name in enumerate(Full_List, start=1):
        file_name = 'Sphere_Results' + str(output_index).zfill(4) + '.vtp'
        output_path = 'STLMerged/' + file_name

        if os.path.isfile(output_path) and os.path.getsize(output_path) > 0:
            file_bounds, file_count = _read_vtp_bounds(output_path)
            bounds = _merge_bounds(bounds, file_bounds)
            number_of_spheres += file_count
        else:
            time_path = 'bodiesInfo/' + time_name
            info_names = sorted(
                name for name in os.listdir(time_path)
                if name.startswith('body') and name.endswith('.info')
            )

            spheres = []
            for info_name in info_names:
                sphere = _read_sphere(time_path + '/' + info_name)
                if sphere is not None:
                    spheres.append(sphere)

            spheres.sort(key=lambda sphere: sphere[0])
            for sphere in spheres:
                bounds = _include_sphere(bounds, sphere)
            number_of_spheres += len(spheres)

            _write_sphere_vtp(output_path, spheres)
            print("-- Reading: ", time_name, flush=True)

        records.append((time_name, file_name))

    if number_of_spheres == 0:
        raise SystemExit(
            'ERROR: no sphere { position; radius; } records were found'
        )

    pvd_path = _write_pvd(records)
    _write_paraview_state(pvd_path, bounds)


if not _recorded_stl_exists():
    _merge_implicit_spheres()
    raise SystemExit(0)


if(not (os.path.isdir('STLMerged'))):
    os.system('mkdir STLMerged')
    time_iter = 0

    for item in Full_List:
        Full_List_II= getParticles_List('bodiesInfo/'+item+'/stlFiles')
        time_iter += 1
        for item_II in Full_List_II:
            with open('bodiesInfo/'+item+'/stlFiles/'+item_II+'.stl', 'r') as file:
                data = file.readlines()    
            with open('STLMerged/STL_Results'+str(time_iter).zfill(4)+'.stl', 'a+') as file:
                file.writelines('solid '+item_II+'.stl'+'\n')
                for i in range(1, len(data) -1):
                    file.writelines(data[i])
                file.writelines('endsolid '+item_II+'\n')
        print("-- Reading: ", item)

elif(os.path.isdir('STLMerged')):
    reduced_list = os.listdir('STLMerged/')
    reduced_list2 = os.listdir('bodiesInfo/')
    for item in reduced_list2:
        if canBeConvertedToFloat(item):
            if(float(item) == int(float(item))):
                reduced_list2[reduced_list2.index(item)] = int(float(item))
            else:
                reduced_list2[reduced_list2.index(item)] = float(item)
    reduced_list2.sort(key = float)
    reduced_list2 = [str(item) for item in reduced_list2]             

    print(" -- Warning: STLMerged/ already exists. The program will append the data to the existing files size: ",len(reduced_list))  
    print(" -- Simulation length: ",len(reduced_list2))  
    time_iter = len(reduced_list)

    print(reduced_list2[len(reduced_list):])

    for item in reduced_list2[len(reduced_list):]:
        Full_List_II= getParticles_List('bodiesInfo/'+item+'/stlFiles')
        time_iter += 1
        for item_II in Full_List_II:
            with open('bodiesInfo/'+item+'/stlFiles/'+item_II+'.stl', 'r') as file:
                data = file.readlines()    
            with open('STLMerged/STL_Results'+str(time_iter).zfill(4)+'.stl', 'a+') as file:
                file.writelines('solid '+item_II+'.stl'+'\n')
                for i in range(1, len(data) -1):
                    file.writelines(data[i])
                file.writelines('endsolid '+item_II+'\n')
        print("-- Reading: ", item)
