# Initial data

Q:=Quiver(["v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13"],
[["v1","v2","a"],["v2","v3","b"],["v3","v4","c"],["v4","v5","d"],["v5","v6","e"],["v6","v7","f"],["v7","v8","g"],["v8","v9","h"],["v9","v10","i"],["v10","v11","j"],["v11","v12","k"],["v12","v13","l"]]);
KQ:=PathAlgebra(GF(47),Q);
AssignGeneratorVariables(KQ);
relns:=[a*b*c*d*e*f*g*h*i,d*e*f*g*h*i*j*k*l];
A:=KQ/relns;
cat := CatOfRightAlgebraModules(A);

zero := ZeroModule(A); 
I1:=IndecInjectiveModules(A)[1];
I2:=IndecInjectiveModules(A)[2];
I3:=IndecInjectiveModules(A)[3];
I4:=IndecInjectiveModules(A)[4];
I5:=IndecInjectiveModules(A)[5];
I6:=IndecInjectiveModules(A)[6];
I7:=IndecInjectiveModules(A)[7];
I8:=IndecInjectiveModules(A)[8];
I9:=IndecInjectiveModules(A)[9];
I10:=IndecInjectiveModules(A)[10];
I11:=IndecInjectiveModules(A)[11];
I12:=IndecInjectiveModules(A)[12];

P1:=IndecProjectiveModules(A)[1];
P2:=IndecProjectiveModules(A)[2];
P3:=IndecProjectiveModules(A)[3];
P4:=IndecProjectiveModules(A)[4];
P5:=IndecProjectiveModules(A)[5];
P6:=IndecProjectiveModules(A)[6];
P7:=IndecProjectiveModules(A)[7];
P8:=IndecProjectiveModules(A)[8];
P9:=IndecProjectiveModules(A)[9];
P10:=IndecProjectiveModules(A)[10];
P11:=IndecProjectiveModules(A)[11];
P12:=IndecProjectiveModules(A)[12];


S1:=SimpleModules(A)[1];
S2:=SimpleModules(A)[2];
S3:=SimpleModules(A)[3];
S4:=SimpleModules(A)[4];
S5:=SimpleModules(A)[5];
S6:=SimpleModules(A)[6];
S7:=SimpleModules(A)[7];
S8:=SimpleModules(A)[8];


M12:=Range(AlmostSplitSequence(S1)[1]);
M23:=Range(AlmostSplitSequence(S2)[1]);
M34:=Range(AlmostSplitSequence(S3)[1]);
M45:=Range(AlmostSplitSequence(S4)[1]);
M56:=Range(AlmostSplitSequence(S5)[1]);
M67:=Range(AlmostSplitSequence(S6)[1]);
M78:=Range(AlmostSplitSequence(S7)[1]);



# Currently not available in QPA 22.05.2025

homologyOfComplex := function( C, i )
    return CoKernel( LiftingInclusionMorphisms( KernelInclusion( DifferentialOfComplex( C, i ) ), ImageInclusion( DifferentialOfComplex( C, i+1 ) ) ) );
end;

goodTruncationBelow := function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist, kerinc, imageproj, imageinc;

    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := PositivePartFrom( difflist, i+2 );
    kerinc := KernelInclusion( DifferentialOfComplex( C, i ) );
    imageproj := ImageProjection( DifferentialOfComplex( C, i+1 ) );
    imageinc := ImageInclusion( DifferentialOfComplex( C, i+1 ) );
    newpart := FiniteInfList( i, [ cat.zeroMap( Source(kerinc), cat.zeroObj ),
                                   imageproj * LiftingInclusionMorphisms( kerinc, imageinc ) ] );
    zeropart := NegativePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i-1 );
    newdifflist := InfConcatenation( truncpart, newpart, zeropart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end;

goodTruncationAbove := function( C, i )
    local cat, difflist, truncpart, newpart, zeropart, newdifflist, factor, factorinclusion;

    cat := CatOfComplex( C );
    difflist := DifferentialsOfComplex( C );
    truncpart := NegativePartFrom( difflist, i-1 );
    factorinclusion := ImageInclusion( DifferentialOfComplex( C, i ) );
    factor := Range( factorinclusion );
    newpart := FiniteInfList( i, [ factorinclusion, cat.zeroMap( cat.zeroObj, factor ) ] );

    zeropart := PositivePartFrom( DifferentialsOfComplex( ZeroComplex( cat ) ),
                                  i+2 );
    newdifflist := InfConcatenation( zeropart, newpart, truncpart );
    
    return ComplexByDifferentialList( cat, newdifflist );

end;

## Functions in project


taum := function(complex, m)
	local tauInverse, projResTauInverse;
	tauInverse := Shift(ProjectiveToInjectiveComplex(BrutalTruncationAbove(ProjectiveResolutionOfComplex(complex), m)), 1);
	projResTauInverse := ProjectiveResolutionOfComplex(tauInverse);
	Print(tauInverse, "\n");
  return projResTauInverse;
end;

taumNew := function(complex, m)
	local tauInverse, projResTauInverse,projRes,projResBound,nakProjResBound,shiftedTemp,tauInv;
	projRes := ProjectiveResolutionOfComplex(complex);
	projResBound := BrutalTruncationAbove(projRes,m);
	nakProjResBound := ProjectiveToInjectiveComplex(projResBound);
	shiftedTemp := Shift(nakProjResBound,1);
	tauInv := goodTruncationBelow(shiftedTemp,0);
	tauInverse := Shift(ProjectiveToInjectiveComplex(BrutalTruncationAbove(ProjectiveResolutionOfComplex(complex), m)), 1);
	projResTauInverse := ProjectiveResolutionOfComplex(tauInverse);
	Print(tauInverse, "\n");
	Print(projResTauInverse, "\n");
  return tauInv;
end;


extendedInjective := function(A, i, m)
  local cat;
  cat := CatOfRightAlgebraModules(A);
  return StalkComplex(cat, IndecInjectiveModules(A)[i], m-1);
end;

#Temporary function for cyclic nakayama m=2
middleOfSequenceTemp := function(x,y)
	local firstComplex, secondComplex, f, chainMap, cone, outputComplex, hom1, hom2, obc0, obc1, obc2,hComp0,hComp1,tempVar;
	
	hom1 := HomOverAlgebra(x[1],x[2]);
	hom2 := HomOverAlgebra(y[1],y[2]);
	
	if x[1]=ZeroModule(A) then
		firstComplex := StalkComplex(cat,x[2],-1);
	elif x[2]=ZeroModule(A) then
		firstComplex := StalkComplex(cat,x[1],0);
	else
		firstComplex := FiniteComplex(cat,0,hom1);
	fi;
	
	if y[1]=ZeroModule(A) then
		secondComplex := StalkComplex(cat,y[2],0);
	elif y[2]=ZeroModule(A) then
		secondComplex := StalkComplex(cat,y[1],1);
	else
		secondComplex := FiniteComplex(cat,1,hom2);
	fi;
	
	f:= HomOverAlgebra(x[1],y[2]);
	chainMap := FiniteChainMap(firstComplex,secondComplex,0,f);
	cone := MappingCone(chainMap);
	outputComplex := cone[1];
	
	tempVar := ObjectOfComplex(outputComplex,2);
	
	#Print(secondComplex);

	obc0 := List(DecomposeModuleViaCharPoly(ObjectOfComplex(outputComplex,0)), a-> DimensionVector(a));
	obc1 := List(DecomposeModuleViaCharPoly(ObjectOfComplex(outputComplex,1)), a-> DimensionVector(a));

	hComp0 := List(DecomposeModuleViaCharPoly(homologyOfComplex(outputComplex,0)), a-> DimensionVector(a));
	hComp1 := List(DecomposeModuleViaCharPoly(homologyOfComplex(outputComplex,1)), a-> DimensionVector(a));
	
	Print("Summands degree 0: ", obc0,"\n");
	Print("Summands degree 1: ",obc1,"\n");
	
	Print("0th homology: ",hComp0,"\n");
	Print("1st homology: ",hComp1,"\n");
	
	return outputComplex;
end;

printFirstHomologies := function(C)
	local h0,h1,h2,h3;
	
	#h0 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,0)), a-> DimensionVector(a));
	#h1 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,1)), a-> DimensionVector(a));
	#h2 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,2)), a-> DimensionVector(a));
	#h3 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,3)), a-> DimensionVector(a));
	
	h0 := DimensionVector(homologyOfComplex(C,0));
	h1 := DimensionVector(homologyOfComplex(C,1));
	h2 := DimensionVector(homologyOfComplex(C,2));
	h3 := DimensionVector(homologyOfComplex(C,3));
	
	Print("0th homology: ",h0,"\n");
	Print("1st homology: ",h1,"\n");
	Print("2nd homology: ",h2,"\n");
	#Print("3rd homology: ",h3,"\n");
	
	return C;
end;