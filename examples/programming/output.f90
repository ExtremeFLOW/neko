! This tutorial demonstrates how to output simulation data to .fld and .csv
! files in Neko. It covers:
! 
! - Writing single fields and lists of fields to .fld files using `fld_file_t`.
! - Controlling output precision for .fld files.
! - Writing data to .csv files using `csv_file_t` and `vector_t`.
! - Managing output across simulation steps and cleaning up resources.
module user
  use neko
  implicit none

  ! A custom field we will output
  type(field_t), target :: my_field1
  type(field_t), target :: my_field2

  ! A custom vector for CSV output
  type(vector_t) :: vec

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%user_startup => user_startup
    user%user_init_modules => user_init_modules
    user%user_finalize_modules => user_finalize_modules

  end subroutine user_setup

  ! We will use the user_startup routine to manipulate the end time.
  subroutine user_startup(params)
    type(json_file), intent(inout) :: params

    call params%add("case.end_time", 0.0_rp)
  end subroutine user_startup

  subroutine user_init_modules(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    ! A writer for .fld files
    type(fld_file_t) :: fld_writer
    ! Storage for a list of fields
    type(field_list_t) :: field_pair
    ! A writer for CSV files
    type(csv_file_t) :: csv_writer

    type(field_t), pointer :: field_ptr1, field_ptr2

    !
    ! WRITING FLD FILES
    ! 

    ! Initialize the fields and set the values to 1.0 and 2.0
    call my_field1%init(coef%dof, "my_field1")
    call my_field2%init(coef%dof, "my_field2")
    call field_cfill(my_field1, 1.0_rp)
    call field_cfill(my_field2, 2.0_rp)

    ! Initialize the writer using the filename
    call fld_writer%init("my_output.fld")

    ! Here is how we can output a single field to the fld file.
    ! The field is written as the pressure field.
    call fld_writer%write(my_field1)

    ! We can also pack more fields inside a single file. To that end we need to
    ! pass a field_list_t object to the writer. The field_list_t is a list of
    ! field pointers.

    ! Initialize the field_list_t object with the number of fields we want to 
    ! pack.
    call field_pair%init(2)
    ! Assign the fields to the list, the first parameter is the index of the
    ! field in the list, the second is the field itself.

    field_ptr1 => my_field1
    field_ptr2 => my_field2 
    call field_pair%assign(1, field_ptr1)
    call field_pair%assign(2, field_ptr2)

    ! Write the list of fields to a different file, they will be put in the 
    ! pressure and temperature inside the fld.
    call fld_writer%init("my_output2.fld")
    call fld_writer%write(field_pair)

    ! If one has 3 fields in the list, they will instead be writtent as velocity
    ! components. Further adding 1 will populate the pressure, then the next
    ! field will be the temperature, and all consecutive fields wil be scalars.

    ! Finally, note that you can set the precision of the output by calling
    call fld_writer%set_precision(dp) ! <- sets double precision output
    call fld_writer%set_precision(sp) ! <- sets single precision output

    !
    ! WRITING CSV FILES
    ! 

    ! A vector_t can be used to pass values to a CSV file writer. Let us init
    ! one to have 5 elements and set them to 3.0.
    call vec%init(5)
    vec = 3.0_rp

    ! We initialize the CSV writer just like the fld writer.
    call csv_writer%init("my_output.csv")
    ! The writer can optionally add a header string to the file.
    call csv_writer%set_header("# p1, p2, p3, p4, p5")

    ! This will write one line to the CSV.
    call csv_writer%write(vec)

    ! Let us add another line. A common scenario is that you will update a
    ! vector at each time step, e.f. in user_check and then append the result
    ! to the CSV file. 
    vec = 4.0_rp
    call csv_writer%write(vec)

    ! Note that the CSV writer will append to file, if it exists, so you have
    ! to take care of that manually. Finally, note that one can also use
    ! matrix_t objects to write to CSV files, in case you have 2D data. 

  end subroutine user_init_modules

  subroutine user_finalize_modules(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! Don't forget to free the types you initialized
    call my_field1%free()
    call my_field2%free()
    call vec%free()
  end subroutine user_finalize_modules

end module user
