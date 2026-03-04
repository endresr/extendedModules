# Initial data

Q:=Quiver(["v1","v2"],
[["v1","v2","a"],["v2","v1","b"]]);

KQ:=PathAlgebra(GF(47),Q);
AssignGeneratorVariables(KQ);
relns:=[a*b,b*a];
A:=KQ/relns;
cat := CatOfRightAlgebraModules(A);

zero := ZeroModule(A); 
I_1:=IndecInjectiveModules(A)[1];
I_2:=IndecInjectiveModules(A)[2];

P_1:=IndecProjectiveModules(A)[1];
P_2:=IndecProjectiveModules(A)[2];

S_1:=SimpleModules(A)[1];
S_2:=SimpleModules(A)[2];

##Showing problem
Print("\n");
Print(TextAttr.bold,"Algebra: Cyclic Nakayama with two nodes and relations of length 2",TextAttr.normal,"\n\n");
Print("Problem can be seen on the simples. ","\n");
Print("Take the projective resolution of simple at 1:","\n\n");
projResS1:=ProjectiveResolutionOfComplex(StalkComplex(cat,S_1,0));
ObjectOfComplex(projResS1,3);
Print(projResS1, "\n\n");
Print("The labels are wrong, but differentials seems to be correct,\n e.g. for differential at second position","\n\n");
Print(TextAttr.bold,"Kernel is:",TextAttr.normal,"\n");
f:=DifferentialOfComplex(projResS1,2);
kerF:=Kernel(f);
Print(kerF,"\n\n");
Print(TextAttr.bold,"Image is:",TextAttr.normal,"\n");
imF:=Image(f);
Print(imF,"\n\n");
Print("However, after applying the Nakayama functor, \n","all differentials are set to identity.","\n\n");
Print("Nakayama functor applied:","\n\n");
nakOfprojResS1:=ProjectiveToInjectiveComplex(projResS1);
Print(nakOfprojResS1,"\n\n");
Print("Look at differential at position 2:","\n\n");
g:=DifferentialOfComplex(nakOfprojResS1,2);
Print(TextAttr.bold,"Kernel is:",TextAttr.normal,"\n");
kerG:=Kernel(g);
Print(kerG,"\n\n");
Print(TextAttr.bold,"Image is:",TextAttr.normal,"\n");
imG:=Image(g);
Print(imG,"\n");


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

printFirstHomologies := function(C)
	local h0,h1,h2,h3,h4;
	
	#h0 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,0)), a-> DimensionVector(a));
	#h1 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,1)), a-> DimensionVector(a));
	#h2 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,2)), a-> DimensionVector(a));
	#h3 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,3)), a-> DimensionVector(a));
	#h3 := List(DecomposeModuleViaCharPoly(homologyOfComplex(C,4)), a-> DimensionVector(a));
	
	h0 := DimensionVector(homologyOfComplex(C,0));
	h1 := DimensionVector(homologyOfComplex(C,1));
	h2 := DimensionVector(homologyOfComplex(C,2));
	h3 := DimensionVector(homologyOfComplex(C,3));
	h4 := DimensionVector(homologyOfComplex(C,4));
	
	Print("0th homology: ",h0,"\n");
	Print("1st homology: ",h1,"\n");
	Print("2nd homology: ",h2,"\n");
	Print("3rd homology: ",h3,"\n");
	Print("4th homology: ",h4,"\n");
	
	return C;
end;

