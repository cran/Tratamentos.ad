#' @importFrom crayon green bold
# --------------------------------------



.onAttach <- function(libname, pkgname) {
  packageStartupMessage(green("############################################################"))
  packageStartupMessage(green(paste0("Obrigado por utilizar o",bold(" Tratamentos.ad "))))
  packageStartupMessage(green("Author: Alcinei Mistico Azevedo (ICA-UFMG)"))
  packageStartupMessage(green("Veja tutoriais sobre este e outros pacotes no youtube:"))
  packageStartupMessage(green("https://www.youtube.com/channel/UCDGyvLCJnv9RtTY1YMBMVNQ"))
  packageStartupMessage(green("Se inscreva e compartilhe para ajudar o canal a crescer."))
  packageStartupMessage(green("############################################################"))

}



