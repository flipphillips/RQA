(* Paclet Info File *)

(* created 2018/03/05*)

Paclet[
    Name -> "RQA",
    Description -> "Recurrance Quantification Package",
    Creator -> "Flip Phillips <flip@skidmore.edu>",
    Publisher -> "Skidmore Vision Lab",
    Copyright -> "(c) 2018- Flip Phillips",
    License -> "MIT",
    Version -> "0.1.1",
    BuildNumber -> "10",
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
        Context -> {"FPTools`"}
      },

      {"FrontEnd", 
        Prepend -> True
      }
    }
]
