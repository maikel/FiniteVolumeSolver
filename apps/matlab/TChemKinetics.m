function [config] = TChemKinetics(mechanism)
  % These are constants used to dispatch certain functions in the .mex file.
  %  kSetMechanism,
  %  kSetPressureIsentropic,
  %  kSetTPX,
  %  kRecPUTY,
  %  kAdvance,
  %  kGetMolarMasses,
  %  kGetSpeciesNames,
  %  kFindSpeciesIndex,
  %  kGetEnthalpies,
  %  kUnknown
  kSetMechanism = 0;
  kSetIsentropicPressure = kSetMechanism + 1;
  kSetTPX = kSetIsentropicPressure + 1;
  kRecSetPUTY = kSetTPX + 1;
  kAdvance = kRecSetPUTY + 1;
  kGetMolarMasses = kAdvance + 1;
  kSpeciesNames = kGetMolarMasses + 1;
  kSpeciesIndex = kSpeciesNames + 1;
  kGetEnthalpies = kSpeciesIndex + 1;
  kSetOdeSolver = kGetEnthalpies + 1;

  TChemKinetics_(kSetMechanism, mechanism);
  
  % load default config
  config = BaseConfig(1);
  config.kineticsif.speciesNamesCellArray = TChemKinetics_(kSpeciesNames);
  
  % Define matlab functions which can be hooked with the library
  
  function [states] = SetTPX(T, P, X)
    states = TChemKinetics_(kSetTPX, T, P, X);
  end
  
  function [states] = SetIsentropicPressure(states, P)
    states = TChemKinetics_(kSetIsentropicPressure, states, P);
  end
  
  function [states, rates] = Advance(states, dt)
    [states, rates] = TChemKinetics_(kAdvance, states, dt);
  end
  
  function [states] = RecSetPUTY(states)
    states = TChemKinetics_(kRecSetPUTY, states);
  end

  function [names] = SpeciesNames()
    names = config.kineticsif.speciesNamesCellArray;
  end

  function [index] = SpeciesIndex(name)
    index = TChemKinetics_(kSpeciesIndex, name) + 1;
  end

  function [enthalpies] = BaseEnthalpies(T)
		if ~exist('T', 'var')
			T = 298.15;
		end
		enthalpies = TChemKinetics_(kGetEnthalpies, T);
  end
     
	function [mw] = GetMolarMasses()
		mw = TChemKinetics_(kGetMolarMasses);
  end

  function SetOdeSolver(solver)
     TChemKinetics_(kSetOdeSolver, solver);
  end

  function SetNumberOfThreads(n_threads) 
    if (n_threads > 1) 
      error('TChem does not support multi-threading.');
    end
  end

  % Set hooks and register the thermo dynamic functions with the library
  config.kineticsif.setTPX = @SetTPX;
  config.kineticsif.recSetPUTY = @RecSetPUTY;
  config.kineticsif.setIsentropicP = @SetIsentropicPressure;
  config.kineticsif.advance = @Advance;
  config.kineticsif.speciesNames = @SpeciesNames;
  config.kineticsif.speciesIndex = @SpeciesIndex;
	config.kineticsif.baseEnthalpies = @BaseEnthalpies;
  config.kineticsif.molarMasses = @GetMolarMasses;
  config.kineticsif.setOdeSolver = @SetOdeSolver;
  config.kineticsif.setNumberOfThreads = @SetNumberOfThreads;
    
  % Set State Variable Indices
  config.v.RHOY = (config.v.RHOY1:(config.v.RHOY1 + ...
                    length(config.kineticsif.speciesNamesCellArray) - 1));
  config.i.CV = 3;
  config.i.CP = 2;
  config.i.T = 1;
  config.i.P = 0;
  config.i.COUNT = 4;   
end

