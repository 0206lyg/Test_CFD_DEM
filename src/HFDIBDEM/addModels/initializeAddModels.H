    // Note: if you want to add a new addModel, do it here

    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

    forAll (bodyNames_,nameI)
    {
        word bodyName(bodyNames_[nameI]);
        const dictionary& bodyDict(HFDIBDEMDict_.subDict(bodyName));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Creating addition model for body: " << bodyName << endl;

        // initialize the addition models
        dictionary addModelDict(bodyDict.subDict("bodyAddition"));
        word howToAdd(word(addModelDict.lookup("addModel")));
        word bodyGeom;
        std::unique_ptr<geomModel> bodyGeomModel;
        // check if the immersedDict_ contains bodyGeom
        if (bodyDict.found("bodyGeom"))
        {
            word input = bodyDict.lookup("bodyGeom");
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for " << bodyName
                << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        // Decouple the logical body/addition-model name from the STL basename.
        // Existing dictionaries remain unchanged when stlBaseName is absent.
        word stlBaseName(bodyName);
        if (bodyDict.found("stlBaseName"))
        {
            stlBaseName = word(bodyDict.lookup("stlBaseName"));
        }
        word stlPath("constant/triSurface/"+stlBaseName+".stl");

        if(bodyGeom == "convex")
        {
            InfoH << iB_Info << "Using STL source for " << bodyName
                << ": " << stlPath << endl;
            bodyGeomModel = std::unique_ptr<convexBody>(
                new convexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "nonConvex")
        {
            InfoH << iB_Info << "Using STL source for " << bodyName
                << ": " << stlPath << endl;
            bodyGeomModel = std::unique_ptr<nonConvexBody>(
                new nonConvexBody(mesh_,stlPath,thrSurf));
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("startPosition");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel = std::unique_ptr<sphereBody>(
                new sphereBody(mesh_,startPosition,radius,thrSurf));
        }
        else
        {
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            InfoH << iB_Info << "Using STL source for " << bodyName
                << ": " << stlPath << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::unique_ptr<nonConvexBody>(
                new nonConvexBody(mesh_,stlPath,thrSurf));
        }

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "trying to set addModel: " << howToAdd << endl;

        if (howToAdd == "once")
        {
            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "onceFromFile")
        {
            addModels_.set
            (
                nameI,
                new addModelOnceFromFile
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_,
                    bodyGeom,
                    thrSurf
                )
            );
        }
        else if (howToAdd == "onceScatter")
        {
            addModels_.set
            (
                nameI,
                new addModelOnceScatter
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "repeatSamePosition")
        {
            addModels_.set
            (
                nameI,
                new addModelRepeatSamePosition
                (
                    addModelDict,
                    mesh_,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "repeatRandomPosition")
        {
            addModels_.set
            (
                nameI,
                new addModelRepeatRandomPosition
                (
                    addModelDict,
                    mesh_,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else if (howToAdd == "distribution")
        {
            addModels_.set
            (
                nameI,
                new addModelDistribution
                (
                    addModelDict,
                    mesh_,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else
        {
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "No addition model specified, defaulting to addModelOnce"
                << endl;

            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
    }
