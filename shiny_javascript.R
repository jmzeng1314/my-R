
## https://github.com/daattali/shinyjs
## 这个包引入的这个toggle函数，可以把某个带ID的div区域给隐藏掉
library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(),
  dashboardBody(
    useShinyjs(),
    actionButton("button", "Click me"),
    div(id = "be_hidden", "Hello!")
  )
)

server <- function(input, output) {
  observeEvent(input$button, {
    toggle("be_hidden")
  })
}

shinyApp(ui, server)
## 下面的代码实现同样的功能：
shinyjs::hidden(
  div(id = "advanced",
    ...
))


## 下面的代码首先给一个带有很多控件的界面
library(shiny)
shinyApp(
  ui = fluidPage(
    div(id = "myapp",
      h2("shinyjs demo"),
      textInput("name", "Name", ""),
      numericInput("age", "Age", 30),
      textInput("company", "Company", ""),
      p("Timestamp: ", span(date())),
      actionButton("submit", "Submit")
    )
  ),

  server = function(input, output) {
  }
)

## 然后通过在server端添加一系列js事件来操作那些控件。
## 第一个是：控制button是否可以被点击
observe({
  ## 只有当用户填写了name的时候，才能按提交这个按钮
  if (is.null(input$name) || input$name == "") {
    shinyjs::disable("submit")
  } else {
    shinyjs::enable("submit")
  }
})

## 下面的代码实现同样的功能
observe({
  shinyjs::toggleState("submit", !is.null(input$name) && input$name != "")
})

## 或者利用下面的代码修改某个ID的CSS属性
checkboxInput("big", "Bigger text", FALSE)
shinyjs::inlineCSS(list(.big = "font-size: 2em"))

## 或者根据控件来调整整个页面的class属性
observe({
  if (input$big) {
    shinyjs::addClass("myapp", "big")
  } else {
    shinyjs::removeClass("myapp", "big")
  }
})
## 与下面代码是同样功能：
observe({
  shinyjs::toggleClass("myapp", "big", input$big)
})

## 或者在console端显示信息
observeEvent(input$submit, {
  shinyjs::info("Thank you!")
})

## 比较复杂而且成熟的例子如下：
library(shiny)
library(shinyjs)

jsCode <- "shinyjs.pageCol = function(params){$('body').css('background', params);}"
## 内部引入一个js语句，如果js很多，需要用文件来引入

shinyApp(
  ui = fluidPage(
    useShinyjs(),
    extendShinyjs(text = jsCode),
    selectInput("col", "Colour:",
                c("white", "yellow", "red", "blue", "purple"))
  ),
  server = function(input, output) {
    observeEvent(input$col, {
      js$pageCol(input$col)
    })
  }
)
## 另外一种引入js代码的方法是：
jscode <- "
shinyjs.init = function() {
  $(document).keypress(function(e) { alert('Key pressed: ' + e.which); });
}"

shinyApp(
  ui = fluidPage(
    useShinyjs(),
    extendShinyjs(text = jscode),
    "Press any key"
  ),
  server = function(input, output) {}
)

如果要引入外部js文件
1、create www folder in the same folder as server.R and ui.R
2、put javascript file into www folder.
3、put tags$head(tags$script(src="hoge.js")) in UI.

或者：tags$head(HTML("<script type='text/javascript' src='js/hoge.js'></script>"))


