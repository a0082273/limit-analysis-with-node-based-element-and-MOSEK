def limt_analysis():
    from openfile import InputFile
    from inputdata import InputData
    import materials
    import mesh
    import libraries

    infile=InputFile()
    infile.open_file()
    infile.get_number_of_line()

    print('')
    print('***************Input file data***************')
    indata=InputData(infile)
    indata.get_title()
    indata.get_problem()
    indata.get_output_file_type()

    mat=materials.Material()
    mat.get_material_type(infile)
    mat.get_material_parameters()

    el=mesh.Element()
    el.get_element_type(infile)
    el.get_number_of_element(infile)
    el.get_node_of_element(infile)

    poi=mesh.Point()
    poi.get_number_of_point(infile)
    poi.get_coordinate(infile)
    poi.get_prescribed_velocity(infile)
    poi.get_prescribed_force(infile)
    poi.get_applied_force(infile)
    print('*********************************************')
    print('')

    primal=libraries.Primal(mesh,materials)
    print('')
    print('***************Primal problem****************')
    primal.assemble()
    primal.solve()
    print('*********************************************')
    print('')

    dual=libraries.Dual(mesh,materials)
    print('')
    print('****************Dual problem*****************')
    dual.assemble()
    dual.solve()
    print('*********************************************')
    print('')

    libraries.output(infile,mesh,mat)

    print('Finished!')
    return

if(__name__=='__main__'):
    limt_analysis()
