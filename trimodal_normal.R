#----------------------------------PACKAGES--------------------------------------
if(!require(hypergeo)){install.packages('hypergeo')} # (fornece a func¸˜ao hipergeom´etrica)
if(!require(pracma)){install.packages('pracma')} # (fornece as duas formas da func¸˜ao gama incompleta)
if(!require(gridExtra)){install.packages('gridExtra')}
if(!require(ggplot2)){install.packages('ggplot2')}

library(hypergeo)
library(pracma)
library(ggplot2)
library(gridExtra)
#-----------------------------ASSISTANT FUNCTIONS--------------------------------
gama.inc <- function(x,p){gamma(p)*pgamma(x,p,1,lower=TRUE)}

# func¸˜ao gama incompleta inferior regularizada, nome "T" para bater com o paper
T <- function(x,alfa,p){ gama.inc(x^2/alfa^2,p) / gamma(p)}

# fator normalizante
z.theta <- function(rho,delta,sigma,alfa,p){
  (rho+delta)*sigma - ((2*delta*sigma*alfa*gamma(p+1/2))/(sqrt(2*pi)*gamma(p))) * hypergeo(1/2,p+1/2,3/2,-(alfa^2)/2)
}

# func¸˜ao I(u;a,b,c)
Int <- function(u,a,b,c){
  integral(f= function(x){erf(a*x)*(x^b)*exp(- c^2 * x^2)}, xmin=0, xmax=u)
}

Int <- Vectorize(Int)

#-----------------------------TRIMODAL DENSITY--------------------------------
# densidade trimodal com g(x)=phi(x) (normal trimodal)
trimodal.normal <- function(x,p,mu,sigma,alfa,rho,delta){
  (1/z.theta(rho,delta,sigma,alfa,p)) * (rho + delta * T( (x-mu)/sigma, alfa, p) ) * dnorm((x-mu)/sigma, mean=0, sd=1)
}

#----------------------CUMULATIVE DISTRIBUTION FUNCTION-------------------------
trimodal.normal.FDA <- function(x,p,mu,sigma,alfa,rho,delta){
  if(x<mu){ as.numeric( sigma/z.theta(rho,delta,sigma,alfa,p)
                        * (rho + delta*T((x-mu)/sigma, alfa, p)) * pnorm((x-mu)/sigma) +
                          delta*sigma/z.theta(rho,delta,sigma,alfa,p) * ( 1/2 - ((alfa*gamma(p+1/2))
                                                                                 / (sqrt(2*pi) * gamma(p))) * hypergeo(1/2,p+1/2,3/2,-(alfa^2)/2) ) -
                          delta*sigma/z.theta(rho,delta,sigma,alfa,p) * ( (1/2)
                                                                          * T((x-mu)/sigma, alfa, p) - (1/gamma(p))
                                                                          * Int(((mu-x)/(alfa*sigma)), alfa/sqrt(2), 2*p-1, 1) ) )
  } else { as.numeric( sigma/z.theta(rho,delta,sigma,alfa,p)
                       * (rho + delta*T((x-mu)/sigma, alfa, p)) * pnorm((x-mu)/sigma) +
                         delta*sigma/z.theta(rho,delta,sigma,alfa,p)
                       * ( 1/2 - ((alfa*gamma(p+1/2)) / (sqrt(2*pi) * gamma(p)))
                           * hypergeo(1/2,p+1/2,3/2,-(alfa^2)/2) ) -
                         delta*sigma/z.theta(rho,delta,sigma,alfa,p)
                       * ( (1/2) * T((x-mu)/sigma, alfa, p) + (1/gamma(p))
                           * Int(((x-mu)/(alfa*sigma)), alfa/sqrt(2), 2*p-1, 1) ) )
  } }

trimodal.normal.FDA = Vectorize(trimodal.normal.FDA)

#--------------------------------CHARTS---------------------------------------
# func¸˜oes para gerar os gr´aficos, em vez de repetir c´odigos
grafico.densidade <- function(p,
                              mu,
                              sigma,
                              delta,
                              rho,
                              alfa,
                              lwd,
                              colour,
                              titulo=titulo,
                              xtitulo=expression(x),
                              ytitulo=expression(f(x,theta)), 
                              subtitle=subtitle){
  df <- data.frame(x=seq(-15,15,0.05),
                   y=as.numeric(trimodal.normal(seq(-15,15,0.05),
                                                p=p,
                                                mu=mu,
                                                sigma=sigma,
                                                delta=delta,
                                                rho=rho,
                                                alfa=alfa)))
  ggplot(df,aes(x,y)) +
    geom_line(lwd=lwd, colour=colour,alpha=0.5) +
    theme_classic() +
    labs(title = titulo,
         subtitle=subtitle,
         x = xtitulo,
         y = ytitulo)
}

grafico.fda <- function(p,mu,sigma,delta,rho,alfa,lwd,colour,titulo=titulo,
xtitulo=expression(x),ytitulo=expression(F(x,theta)), subtitle=subtitle){
  df <- data.frame(x=seq(-15,15,0.05),
        y=as.numeric(trimodal.normal.FDA(seq(-15,15,0.05), p=p,mu=mu,
        sigma=sigma, delta=delta, rho=rho, alfa=alfa)))
  
  ggplot(df,aes(x,y)) +
    geom_line(lwd=lwd, colour=colour,alpha=0.5) +
    theme_classic() +
    labs(title = titulo,
         subtitle=subtitle,
         x = xtitulo,
         y = ytitulo)
}

#-------------------------------TRIMODAIS---------------------------------------
g1.densidade = grafico.densidade(p=6,
                  mu=0,
                  sigma=3,
                  delta=4,
                  rho=1.2,
                  alfa=0.6,
                  lwd=1,
                  colour='gray',
                  titulo='Trimodal Normal Density',
                  subtitle=expression(
                    paste(p,'=6, ',mu,'=0, ',sigma,'=3, ', delta,'=4, ',rho,'=1.2, ',alpha,'=0.6')))  

g1.fda <- grafico.fda(p=6,mu=0,sigma=3,delta=4,rho=1.2,alfa=0.6,lwd=1,colour='gray',
                      titulo='Trimodal Normal Distribution Function',
                      subtitle=expression(paste(p,'=6, ',mu,'=0, ',sigma,
                                                '=3, ',delta,'=4, ',rho,'=1.2, ',alpha,'=0.6')))

trimodal = grid.arrange(g1.densidade, g1.fda, nrow=1, ncol=2)

#ggsave('trimodal.pdf', plot=trimodal4, device = cairo_pdf, 
#       path = 'C:\\Users\\04291229186\\Desktop', dpi = 'retina', width=30, height=18, units='cm' )

