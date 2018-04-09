#
# Be sure to run `pod lib lint VitalSignsHandler.podspec' to ensure this is a
# valid spec before submitting.
#
# Any lines starting with a # are optional, but their use is encouraged
# To learn more about a Podspec see http://guides.cocoapods.org/syntax/podspec.html
#

Pod::Spec.new do |s|
  s.name             = 'VitalSignsHandler'
  s.version          = '0.1.1'
  s.summary          = 'A short description of VitalSignsHandler.'

# This description is used to generate tags and improve search results.
#   * Think: What does it do? Why did you write it? What is the focus?
#   * Try to keep it short, snappy and to the point.
#   * Write the description between the DESC delimiters below.
#   * Finally, don't worry about the indent, CocoaPods strips it!

  s.description      = <<-DESC
TODO: Add long description of the pod here.
                       DESC

  s.homepage         = 'https://github.com/MariMiMari/VitalSignsHandler'
  # s.screenshots     = 'www.example.com/screenshots_1', 'www.example.com/screenshots_2'
  s.license          = { :type => 'MIT', :file => 'LICENSE' }
  s.author           = { 'Мария Тимофеева' => 'mi.maritim@gmail.com' }
  s.source           = { :git => 'https://github.com/MariMiMari/VitalSignsHandler.git', :tag => s.version.to_s }
  # s.social_media_url = 'https://twitter.com/<TWITTER_USERNAME>'

  s.ios.deployment_target = '8.0'

  s.libraries =  'c++'
  s.source_files = 'VitalSignsHandler/Classes/**/**'
  # s.vendored_frameworks = 'VitalSignsHandler/VitalSignsHandler.framework' 
  # s.resource_bundles = {
  #   'VitalSignsHandler' => ['VitalSignsHandler/Assets/*.png']
  # }
  s.requires_arc = true
  s.xcconfig = {
     'CLANG_CXX_LANGUAGE_STANDARD' => 'c++11',
     'CLANG_CXX_LIBRARY' => 'libc++' ,
     'SWIFT_VERSION' => '4.0'
  }

  # s.public_header_files = 'Pod/Classes/**/*.h'
  # s.frameworks = 'UIKit', 'MapKit'
  # s.dependency 'AFNetworking', '~> 2.3'
end
