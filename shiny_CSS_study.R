参考：http://shiny.rstudio.com/articles/css.html

两种方式：
一是直接引用css文件！
shinyUI(fluidPage(

  includeCSS("styles.css"), ## 可以是网页url地址，或者本地文件，相对路径。
    
  headerPanel("New Application"),
  
  sidebarPanel(
    sliderInput("obs", "Number of observations:", 
                min = 1, max = 1000, value = 500)
  ),
  
  mainPanel(plotOutput("distPlot"))
))
而是直接用tags$style

shinyUI(fluidPage(

  tags$head(
    tags$style(HTML("
      @import url('//fonts.googleapis.com/css?family=Lobster|Cabin:400,700');
      
      h1 {
        font-family: 'Lobster', cursive;
        font-weight: 500;
        line-height: 1.1;
        color: #48ca3b;
      }

    ")) ## end for tags
  ), ## end for head
    
  headerPanel("New Application"),
  
  sidebarPanel(
    sliderInput("obs", "Number of observations:", 
                min = 1, max = 1000, value = 500)
  ),
  
  mainPanel(plotOutput("distPlot"))
))
Add CSS to a Shiny UI just as you would add CSS to a web page.

To get CSS into your Shiny App, you can:

Link to an external CSS file
Include CSS in the header of the web page that the app is built on, or
Pass style information to individual HTML elements in your app.
My examples give a glimpse into the options CSS offers your Shiny App. Explore more at Bootswatch and see what you can create.

## 一个成功的例子
在UI端，我加入下面这个代码：
  tags$style(HTML("
            .bigger_image:hover{
                   width : 358px;
                  height : 358px
                  }
                  
                  ") 
  )
# 一个很简单的css,就是对所有的class等于bigger_image的标签的图片都添加一个动作，就是鼠标放在图片上面，图片就变大了~
genes_img=paste0("<img src='http://~~~~~~~~~~~~'  alt='No boxplot!'   class='bigger_image'  width=100 height=100  ></img>")
##那么这一系列图片都是可以调整大小的啦！




