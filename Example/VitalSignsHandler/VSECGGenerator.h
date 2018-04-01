//
//  SensorGenerator.h
//  MirCode
//
//  Created by Мария Тимофеева on 20.11.2017.
//  Copyright © 2017 ___matim___. All rights reserved.
//

@interface VSECGGenerator : NSObject
@property (strong, nonatomic) NSNumber *value;
-(void)start;
-(void)stop;
@end
