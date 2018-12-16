Subroutine Randomize_Model(Params, Nodes, ModelIn, ModelOut)
  Use Param_structure
  Use Nodes_information
  Use Model_structure
  Implicit None
  Type (Nodes_info) :: Nodes
  Type (Parameters) :: Params
  Type (Model) :: ModelIn, ModelOut
  Integer :: i, j, idx
  Real :: kk
  Logical, save :: firsttime = .TRUE.
!
  If(firsttime) Then
     Call RANDOM_SEED
     firsttime = .FALSE.
  End if
!
  ModelOut=ModelIn

  If (Nodes%n_nodes_t .gt. 0) then
     Call Random_number(kk)
     ModelOut%temp=ModelOut%temp + (kk-0.5)*600.
  End if

  If (Nodes%n_nodes_v .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_los=ModelOut%v_los + (kk-.5)*5.e5
  End if

  If (Nodes%n_nodes_mic .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_mic=ModelOut%v_mic + kk*3.e5
  End if

  If (Nodes%n_nodes_blong .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_long=ModelOut%b_long + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_bx .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_x=ModelOut%b_x + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_by .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_y=ModelOut%b_y + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_stray .gt. 0) then
     Call Random_number(kk)
     ModelOut%stray=kk*.7
  End if

  If (Nodes%n_nodes_ffactor .gt. 0) then
     Call Random_number(kk)
     ModelOut%ffactor=kk*.5
  End if

  If (Nodes%n_nodes_chrom_x .gt. 0) then
     Call Random_number(kk)
     ModelOut%ffactor=-5+(kk-.5)
  End if

  If (Nodes%n_nodes_chrom_y .gt. 0) then
     Call Random_number(kk)
     ModelOut%ffactor=1000.+(kk-.5)*500
  End if

  If (Nodes%n_nodes_spic_temp .gt. 0) then
     Call Random_number(kk)
     ModelOut%spic_temp=2000.+(kk)*5000
  End if

  If (Nodes%n_nodes_spic_dens_factor .gt. 0) then
     Call Random_number(kk)
     ModelOut%spic_temp=.1+kk*2
  End if

  If (Nodes%n_nodes_mac .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_mac=kk*2.e5
  End if

  If (Nodes%n_nodes_ab .gt. 0) then
     Do idx=1, Nodes%n_nodes_ab
        Call Random_number(kk)
        ModelOut%abundance(Nodes%i_nodes_ab(idx))= &
             ModelOut%abundance(Nodes%i_nodes_ab(idx))+(kk-0.5)*.2
     End do
  End if

  ! debug
  !  ModelOut%v_mac=ModelIn%v_mac
  ! ModelOut%v_mac=2.e5  
  
End Subroutine Randomize_Model


