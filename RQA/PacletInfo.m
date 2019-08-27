(* Paclet Info File *)

Paclet[
    Name -> "RQA",
    Description -> "Recurrance Quantification Analysis Package",
    Creator -> "Flip Phillips <flip@skidmore.edu>",
    Publisher -> "Skidmore Vision Lab",
    Copyright -> "(c) 2015- Flip Phillips",
    License -> "MIT",
    Version -> "0.3.0",
    BuildNumber -> "26",
    MathematicaVersion -> "10.0+",
    URL -> "https://github.com/flipphillips/RQA",
    Thumbnail -> "Documentation/icon.png",
    (* Loading -> Automatic,     *)
    
    Extensions -> {
      { "Documentation", 
        MainPage -> "Guides/RQA", 
        Language -> "English"},

      { "Kernel", 
        Root -> "Kernel", 
        Context -> {"RQA`"}
      },

      {"FrontEnd", 
        Prepend -> True
      }
    }
]